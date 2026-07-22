#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path


def suppress_first_run_prompt():
    # On a machine that has never run Streamlit before, `streamlit run`
    # blocks on an interactive "enter your email" onboarding prompt before
    # the server starts - it looks like the command just hangs. Streamlit
    # skips that prompt once ~/.streamlit/credentials.toml exists, which is
    # its own documented way to make `streamlit run` non-interactive.
    credentials_path = Path.home() / ".streamlit" / "credentials.toml"
    if credentials_path.is_file():
        return
    credentials_path.parent.mkdir(parents=True, exist_ok=True)
    credentials_path.write_text("[general]\nemail = \"\"\n")


def gui(args):
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print(
            "Streamlit is not installed. Install it with: pip install streamlit",
            file=sys.stderr,
        )
        sys.exit(1)

    suppress_first_run_prompt()

    app_path = Path(__file__).resolve().parent / "streamlit_app.py"
    command = [sys.executable, "-m", "streamlit", "run", str(app_path)]
    if args.port:
        command.extend(["--server.port", str(args.port)])

    subprocess.run(command)
