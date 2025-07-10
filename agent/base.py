# base.py
# -------------------------------

import os
import re
import uuid
from dotenv import load_dotenv
from typing import Annotated, List
from typing_extensions import TypedDict
from langchain_openai import AzureChatOpenAI
from langchain_core.messages import HumanMessage, AIMessage
from langchain_core.tools import tool
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langgraph.graph import StateGraph
from langgraph.graph.message import add_messages
from langgraph.prebuilt import ToolNode, tools_condition
from langgraph.checkpoint.memory import MemorySaver
from datetime import datetime
import traceback
import subprocess
import threading
import time
import json
import tempfile
import pprint
import gradio as gr
from toolkit import tools
from utils import gr_css, get_llm_from_env, chatbot_theme, base_path

# -------------------------------
# Store transcript outside the session
terminal_transcript = ""

def create_cellatria(env_path):

    env_file = os.path.join(env_path, ".env")
    if not os.path.isfile(env_file):
        raise FileNotFoundError(f"*** 🚨 .env file not found at: {env_file}")
    llm = get_llm_from_env(env_path)

    # -------------------------------
    # Define prompt template
    # Load the system message from file

    with open(os.path.join(base_path, "system_prompts.md"), "r") as f:
        system_message = f.read()

    # Defines the agent’s role and toolset
    prompt = ChatPromptTemplate.from_messages([
        ("system", system_message),
        MessagesPlaceholder("messages")
    ])

    # -------------------------------
    # Bind tools to model
    llm_with_tools = llm.bind_tools(tools)
    chat_fn = prompt | llm_with_tools

    # -------------------------------
    # LangGraph schema
    class AgentState(TypedDict):
        messages: Annotated[List, add_messages]

    # -------------------------------
    # Create the graph with schema
    graph_builder = StateGraph(AgentState)
    graph_builder.add_node("tools", ToolNode(tools))
    graph_builder.add_node("chatbot", lambda state: {
        "messages": chat_fn.invoke(state["messages"])
    })
    graph_builder.add_edge("tools", "chatbot")
    graph_builder.add_conditional_edges("chatbot", tools_condition)
    graph_builder.set_entry_point("chatbot")
    graph = graph_builder.compile(checkpointer=MemorySaver())

    # -------------------------------
    # Configuration
    LOG_PATH = "/tmp/cellatria_log.txt"

    # -------------------------------
    # Logging
    def log_status(message: str):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(LOG_PATH, "a") as f:
            # f.write(f"[{timestamp}] {message}\n")
            f.write(f"{message}\n")

    def read_log():
        try:
            with open(LOG_PATH, "r") as f:
                return f.read()
        except FileNotFoundError:
            return "📁 No logs yet."

    # -------------------------------
    # Terminal Session
    class TerminalSession:
        def __init__(self):
            self.process = subprocess.Popen(
                ["bash"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )
            self.lock = threading.Lock()
            self.latest_output = ""
            self.history = ""  

            # Start thread to continuously read output
            threading.Thread(target=self._read_output, daemon=True).start()

        def _read_output(self):
            for line in self.process.stdout:
                with self.lock:
                    self.latest_output += line

        def run_command(self, cmd):
            if self.process.poll() is not None:
                return "❌ Shell has exited."

            with self.lock:
                self.latest_output = ""

            self.process.stdin.write(cmd + "\n")
            self.process.stdin.flush()
            
            time.sleep(0.3)
            
            with self.lock:
                result = self.latest_output.strip()
            return result

    # -------------------------------
    # Initialize the terminal session
    terminal = TerminalSession()

    # -------------------------------
    # Strip ANSI codes before displaying
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

    def clean_ansi(text):
        return ansi_escape.sub('', text)

    # -------------------------------
    # Define the function used in Gradio UI
    def terminal_interface(cmd):
        global terminal_transcript
        raw_output = terminal.run_command(cmd)
        cleaned_output = clean_ansi(raw_output).strip()

        silent_commands = [
            "cd ",       # change directory
            "export ",   # set environment variable
            "unset ",    # unset environment variable
            "alias ",    # define an alias
            "unalias ",  # remove an alias
            "set ",      # set shell options or positional parameters
            "shopt ",    # shell options (bash-specific)
            "trap ",     # set signal handlers
            "pushd ",    # push directory onto stack
            "popd",      # pop directory from stack
            "exec ",     # replace shell with command
            "true",      # does nothing, returns 0
            "false",     # does nothing, returns 1
            ":",         # no-op command
            "clear",     # clears the screen
            "reset",     # resets the terminal
            "wait",      # waits for background jobs
            "disown",    # removes jobs from shell's job table
            "bg",        # resume job in background
            "fg",        # resume job in foreground
            "jobs",      # list background jobs
            "readonly ", # mark variables as read-only
            "declare ",  # variable declaration (bash)
            "typeset ",  # synonym for declare (ksh/zsh)
            "let ",      # arithmetic evaluation
            "source ",   # source a script
            ".",         # shorthand for source
        ]
        is_silent = any(cmd.strip().startswith(sc) for sc in silent_commands)
        
        if cleaned_output:
            terminal_transcript += f"\n$ {cmd}\n{cleaned_output}\n---"
        elif is_silent:
            terminal_transcript += f"\n$ {cmd}\n---"  # Silent: do not show fake "(no output)"
        else:
            # For commands like `ls` or `cat` returning true empty output
            terminal_transcript += f"\n$ {cmd}\n---"

        return terminal_transcript.strip(), ""  # keep input cleared

    # -------------------------------
    # Chat Handler
    # Persistent thread ID for LangGraph
    chat_thread_id = str(uuid.uuid4())
    def gr_block_fn(user_input, pdf_file, history):
        messages = []

        # Ensure initial message is included only once
        if not history:
            history = [initial_message]

        log_status("🟢 New interaction started.")

        # Convert history into LangChain format
        for h in history:
            role = h["role"]
            content = h["content"]
            if role == "user":
                messages.append(HumanMessage(content=content))
            elif role == "assistant":
                messages.append(AIMessage(content=content))

        # Append new user input
        if user_input:
            log_status(f"👤 User input: {user_input}")
            messages.append(HumanMessage(content=user_input))

        # Prepare config
        config = {"configurable": {"thread_id": chat_thread_id}}

        # Invoke LangGraph
        try:
            log_status("🤖 Invoking agent...")
            result = graph.invoke({"messages": messages}, config=config)
            final_message = result["messages"][-1]
            log_status("✅ Agent response received.")
        except Exception as e:
            log_status(f"❌ Error: {str(e)}")
            log_status(traceback.format_exc())
            final_message = AIMessage(content="There was an error processing your request.")

        # If PDF uploaded, include that info too (future use)
        if pdf_file:
            pdf_note = f"\n\n📄 Received PDF: `{pdf_file.name}`. \nI can extract metadata from it!"
            log_status("🟣 Interaction complete.\n---")
            return (
                history + [
                    {"role": "user", "content": user_input},
                    {"role": "assistant", "content": pdf_note}],
                "",
                None,
                history + [
                    {"role": "user", "content": user_input},
                    {"role": "assistant", "content": pdf_note}]
            )
        else:
            pdf_note = ""

        log_status("🟣 Interaction complete.\n---")
        
        return (
            history + [
                {"role": "user", "content": user_input},
                {"role": "assistant", "content": final_message.content + pdf_note}
            ],
            "",
            None,
            history + [
                {"role": "user", "content": user_input},
                {"role": "assistant", "content": final_message.content + pdf_note}
            ]
        )

    # -------------------------------
    # Backend Handler
    def gr_block_stream_backend(user_input, pdf_file, history):
        messages = []
        if not history:
            history = [initial_message]

        for h in history:
            role = h["role"]
            content = h["content"]
            if role == "user":
                messages.append(HumanMessage(content=content))
            elif role == "assistant":
                messages.append(AIMessage(content=content))

        if user_input:
            messages.append(HumanMessage(content=user_input))

        config = {"configurable": {"thread_id": chat_thread_id}}
        backend_log = []

        try:
            for step in graph.stream({"messages": messages}, config=config):
                # --- Extract summary info ---
                summary_lines = []
                if "chatbot" in step:
                    msg = step["chatbot"]["messages"]
                    meta = getattr(msg, "response_metadata", {})
                    model = meta.get("model_name", "")
                    tokens = meta.get("token_usage", {})
                    summary_lines.append(f"**Step:** `chatbot`")
                    summary_lines.append(f"**Model:** `{model}`")
                    summary_lines.append(
                        f"**Tokens:** {tokens.get('completion_tokens', '?')} completion, "
                        f"{tokens.get('prompt_tokens', '?')} prompt, "
                        f"{tokens.get('total_tokens', '?')} total"
                    )
                elif "tools" in step:
                    for tool_msg in step["tools"]["messages"]:
                        summary_lines.append(f"**Step:** `tools`")
                        summary_lines.append(f"**Tool:** `{getattr(tool_msg, 'name', 'unknown')}`")
                # --- End summary info ---

                # Pretty-print the step dict as plain text
                step_txt = pprint.pformat(step, compact=True, width=120)
                backend_log.append(
                    "\n".join(summary_lines) +
                    "\n\n```python\n" + step_txt + "\n```\n"
                )
                yield "\n\n".join(backend_log)
        except Exception as e:
            backend_log.append(f"❌ Error: {str(e)}")
            backend_log.append(traceback.format_exc())
            yield "\n\n".join(backend_log)

    # -------------------------------
    # File Browser Handlers
    def fb_list_subdirs_and_files(path):
        try:
            items = os.listdir(path)
            dirs = sorted([item for item in items if os.path.isdir(os.path.join(path, item))])
            files = sorted([item for item in items if os.path.isfile(os.path.join(path, item))])
            return dirs, files, None
        except Exception as e:
            return [], [], f"❌ Error: {str(e)}"

    def fb_get_dropdown_choices(path):
        dirs, _, _ = fb_list_subdirs_and_files(path)
        dirs = [f"📁 {d}" for d in dirs]
        parent = os.path.dirname(path.rstrip("/"))
        if parent and os.path.abspath(parent) != os.path.abspath(path):
            return [".. (Up)"] + dirs
        return dirs

    def fb_initial_refresh(path):
        choices = fb_get_dropdown_choices(path)
        dirs, files, error = fb_list_subdirs_and_files(path)
        if error:
            file_display_val = f"<span style='color:red'>{error}</span>"
        else:
            file_display_val = "\n".join(f"📄 {f}" for f in files) or "No files in this directory."
        current_path = f"**Current Path:** `{path}`  \n**Folders:** {len(dirs)} | **Files:** {len(files)}"
        return (
            gr.Dropdown(choices=choices, value=None, interactive=bool(choices)),
            file_display_val,
            path,
            current_path
        )

    def fb_navigate_subdir(subdir, base):
        if subdir and subdir.startswith("📁 "):
            subdir = subdir[2:]
        if subdir == ".. (Up)":
            new_path = os.path.dirname(base.rstrip("/"))
            if not new_path:
                new_path = base
        elif subdir:
            new_path = os.path.join(base, subdir)
        else:
            new_path = base
        choices = fb_get_dropdown_choices(new_path)
        dirs, files, error = fb_list_subdirs_and_files(new_path)
        if error:
            file_display_val = f"<span style='color:red'>{error}</span>"
        else:
            file_display_val = "\n".join(f"📄 {f}" for f in files) or "No files in this directory."
        current_path = f"**Current Path:** `{new_path}`  \n**Folders:** {len(dirs)} | **Files:** {len(files)}"
        return (
            gr.Dropdown(choices=choices, value=None, interactive=bool(choices)),
            file_display_val,
            new_path,
            current_path
        )
               
    # -------------------------------
    # Chat Export Feature
    def export_chat_history(state_data):
        uid = uuid.uuid4().hex[:8]
        filename = f"chat_{uid}.json"
        filepath = os.path.join(tempfile.gettempdir(), filename)
        with open(filepath, "w") as f:
            json.dump(state_data, f, indent=2)
        return filepath

    # -------------------------------
    # Initial message
    initial_message = {
        "role": "assistant",
        "content": (
            "👋 Hello! I'm **cellAtria**, your assistant for analyzing single-cell RNA-seq (scRNA-seq) datasets.\n\n"
            "Here's what I can help you with:\n"
            "1. Extract structured metadata from scientific articles (PDF or URL).\n"
            "2. Store and organize metadata in structured project directories.\n"
            "3. Access public databases (currently support GEO) and fetch associated sample metadata.\n"
            "4. Download and organize scRNA-seq datasets.\n"        
            "5. Trigger CellExpress standardized single-cell data processing.\n\n"
            "To see a list of all available actions, type `'help'`.\n"
            "**Let's see how far we can fly together.** 🕊️\n\n"
            "To get started, let's set up your working directory. I can show you the current one or help create a new path. You can change your working directory anytime you wish.\n"
            f"📂 **Current Working Directory:** `{os.getcwd()}`"
        )
    }

    # -------------------------------
    # Clear the log file when app starts
    open("/tmp/cellatria_log.txt", "w").close()  # clears the log file

    # -------------------------------
    # Interface
    with gr.Blocks(theme=chatbot_theme, css=gr_css) as cellatria:
        gr.HTML("""
            <div style='text-align: center; margin-bottom: 0;'>
                <h1 style='margin-bottom: 0.2em;'>Welcome to cellAtria</h1>
                <h3 style='margin-top: 0.2em;'>Agentic Triage of Regulated single cell data Ingestion and Analysis</h3>
            </div>
        """)
        gr.Image(
            value=os.path.join(base_path, "cellatria_header.png"),
            label="Welcome to cellAtria",
            interactive=False,
            show_label=False,
            show_download_button=False,
            container=True,
            show_fullscreen_button=False,
            height=125, 
        )

        # Chat Interface
        chatbot = gr.Chatbot(
            value=[initial_message],
            label="cellAtria Agent",
            show_label=False,
            type="messages",
            height=500,
            show_copy_button=False,
            autoscroll=True,
            resizable=True
        )

        with gr.Row(equal_height=True, elem_id="fixed_top_row"):
            with gr.Column(scale=4):
                user_input = gr.Textbox(placeholder="Give me a task...", label="Your Prompt", lines=1)
            with gr.Column(scale=1):
                with gr.Row():  
                    pdf_upload = gr.File(file_types=[".pdf"], label=".pdf", show_label=True, interactive=True)
                    submit_btn = gr.Button("Submit Prompt/PDF", variant="primary")
            with gr.Column(scale=2):
                log_viewer = gr.Textbox(label="Live Logs", lines=12, interactive=False, elem_id="log_viewer_aes")


        # Hidden state to maintain chat memory
        state = gr.State([initial_message])
        
        # Bind inputs to gr_block_fn
        user_input.submit(
            fn=gr_block_fn,
            inputs=[user_input, pdf_upload, state],
            outputs=[chatbot, user_input, pdf_upload, state]
        )
        submit_btn.click(
            fn=gr_block_fn,
            inputs=[user_input, pdf_upload, state],
            outputs=[chatbot, user_input, pdf_upload, state]
        )

        # --- Backend Panel ---
        with gr.Accordion("Agent Backend", open=False, elem_id="agent_backend_panel"):
            agent_backend_md = gr.Markdown("No agent activity yet.", elem_id="agent_backend_md")

        user_input.submit(
            fn=gr_block_stream_backend,
            inputs=[user_input, pdf_upload, state],
            outputs=agent_backend_md
        )
        submit_btn.click(
            fn=gr_block_stream_backend,
            inputs=[user_input, pdf_upload, state],
            outputs=agent_backend_md
        )

        # --- Terminal Panel ---
        with gr.Accordion("Terminal Panel", open=False, elem_id="logs_terminal_panel"):
            with gr.Row(equal_height=True):
                shell_input = gr.Textbox(placeholder="Enter shell command (e.g., ls -la)", lines=1, label="Command")
            shell_output = gr.Textbox(label="Terminal Output", lines=10, interactive=False, elem_id="terminal_aes")

        # Timer Hook
        log_timer = gr.Timer(value=1.0, active=True)
        log_timer.tick(fn=read_log, inputs=[], outputs=[log_viewer] )

        # Bind inputs to gr_block_fn
        shell_input.submit(
            fn=terminal_interface,
            inputs=shell_input,
            outputs=[shell_output, shell_input]
        )

        # --- File Browser ---
        with gr.Accordion("File Browser", open=False, elem_id="logs_browser_panel"):

            fb_base_path = gr.Textbox(value=os.getcwd(), label="Directory Path")
            fb_current_path_md = gr.Markdown()
            fb_dir_dropdown = gr.Dropdown(label="Subdirectories", choices=[], interactive=True)
            fb_refresh_button = gr.Button("Refresh Directory", variant="primary")
            fb_file_display = gr.Textbox(label="Files in Directory", lines=10, interactive=False)

            fb_dir_dropdown.change(
                fn=fb_navigate_subdir,
                inputs=[fb_dir_dropdown, fb_base_path],
                outputs=[fb_dir_dropdown, fb_file_display, fb_base_path, fb_current_path_md]
            )

            fb_refresh_button.click(
                fn=fb_initial_refresh,
                inputs=fb_base_path,
                outputs=[fb_dir_dropdown, fb_file_display, fb_base_path, fb_current_path_md]
            )

        # --- History Panel ---
        with gr.Accordion("Export Chat", open=False, elem_id="logs_history_panel"):
            export_btn = gr.Button("Download Chat", variant="primary")
            chat_file = gr.File(file_types=[".json"], label=".json", show_label=True, interactive=False)

        # Bind inputs to gr_block_fn
        export_btn.click(
            fn=export_chat_history,
            inputs=[state],
            outputs=[chat_file]
        )   

    return graph, cellatria

# -------------------------------