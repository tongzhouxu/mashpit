import os
import threading
import subprocess
import shutil
import pandas as pd
from flask import Flask, render_template, request, redirect, url_for

app = Flask(
    __name__,
    template_folder=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "templates"
    ),
)

# Lock for thread-safe operations on global variables
lock = threading.Lock()

# Global variables to store the output of the classification and status
mashpit_out = []
mashpit_status = "idle"


def run_mashpit(asm_file_path):
    global mashpit_out
    global mashpit_status
    try:
        with lock:
            mashpit_status = "in progress"
        upload_folder = os.path.join(app.config["UPLOAD_FOLDER"])
        subprocess.run(["mashpit", "query", asm_file_path, upload_folder], check=True)
        query_name = os.path.basename(asm_file_path).split(".")[0]
        # move the output to the static folder
        static_folder = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "static"
        )
        # check if file exists and remove it
        if os.path.exists(os.path.join(static_folder, f"{query_name}_output.csv")):
            os.remove(os.path.join(static_folder, f"{query_name}_output.csv"))
        if os.path.exists(os.path.join(static_folder, f"{query_name}_tree.png")):
            os.remove(os.path.join(static_folder, f"{query_name}_tree.png"))
        shutil.move(f"{query_name}_output.csv", static_folder)
        shutil.move(f"{query_name}_tree.png", static_folder)

        df_output = pd.read_csv(os.path.join(static_folder, f"{query_name}_output.csv"))
        # remove first column
        df_output = df_output.iloc[:, 1:]
        # put column 'asm_acc', 'similarity_score' and 'SNP_tree_link' to the front
        df_output = df_output[
            ["asm_acc", "similarity_score", "SNP_tree_link"]
            + [
                col
                for col in df_output.columns
                if col not in ["asm_acc", "similarity_score", "SNP_tree_link"]
            ]
        ]

        image_path = f"{query_name}_tree.png"
        df_output_html = df_output.to_html(
            index=False, classes="table table-striped table-hover"
        )

        with lock:
            mashpit_out = [image_path, df_output_html]
            mashpit_status = "complete"

    except Exception as e:
        with lock:
            mashpit_status = f"error: {str(e)}"


@app.route("/", methods=["GET"])
def upload():
    return render_template("upload.html")


@app.route("/submit", methods=["POST"])
def submit():
    uploaded_files = request.files
    paths = []
    for key in ["formFile_db", "formFile_sig", "formFile_assembly"]:
        file = uploaded_files[key]
        file_path = os.path.join(app.config["UPLOAD_FOLDER"], file.filename)
        file.save(file_path)
        paths.append(file_path)

    global mashpit_thread
    global mashpit_status

    with lock:
        mashpit_status = "preparing"
        mashpit_thread = threading.Thread(target=run_mashpit, args=(paths[2],))
        mashpit_thread.start()

    return redirect(url_for("loading"))


@app.route("/loading")
def loading():
    global mashpit_status
    with lock:
        current_status = mashpit_status
    return render_template("loading.html", status=current_status)


@app.route("/result")
def result():
    global mashpit_out
    with lock:
        image_path, df_output_html = mashpit_out
    return render_template(
        "result.html", image_path=image_path, df_output_html=df_output_html
    )


@app.route("/status")
def status():
    global mashpit_status
    with lock:
        current_status = mashpit_status
    return current_status


def webserver(args):
    app.config["UPLOAD_FOLDER"] = "uploads"
    # if the folder exists, clear it
    if os.path.exists(app.config["UPLOAD_FOLDER"]):
        for file in os.listdir(app.config["UPLOAD_FOLDER"]):
            os.remove(os.path.join(app.config["UPLOAD_FOLDER"], file))
    else:
        os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)
    app.run(debug=True, host="0.0.0.0", port=8080)
