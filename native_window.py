import webview
import sys

url = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:8501"
webview.create_window(
    "BIOBIO — Molecular Docking & MD Suite",
    url,
    width=1440,
    height=920,
    min_size=(1024, 700),
)
webview.start()
