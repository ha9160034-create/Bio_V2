"""
Bio-Impact Analyzer - Application Entrypoint
============================================
Streamlit entrypoint for standard deployments and Streamlit Community Cloud.
Delegates to the main application in b.py.
"""

from b import main

if __name__ == "__main__":
    main()
