"""
Automated GitHub Uploader for Bio-Impact Analyzer
================================================
Uploads the complete codebase directly to a user's GitHub account via the GitHub REST API.
No local Git installation or administrator privileges required!
"""

import os
import sys
import json
import base64
import argparse
from pathlib import Path

try:
    import requests
except ImportError:
    print("Error: 'requests' library is required. Run via .venv/Scripts/python.")
    sys.exit(1)

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Directories and files to exclude from upload
EXCLUDE_DIRS = {
    ".venv",
    "__pycache__",
    ".git",
    ".pytest_cache",
    ".vscode",
    "md_runs",  # Local simulation outputs
}

EXCLUDE_EXTENSIONS = {
    ".pyc",
    ".pyo",
    ".tmp",
}


def get_files_to_upload(root_dir: Path):
    file_list = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Exclude directories in-place
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS and not d.startswith(".venv")]
        for f in filenames:
            ext = os.path.splitext(f)[1].lower()
            if ext in EXCLUDE_EXTENSIONS:
                continue
            full_path = Path(dirpath) / f
            rel_path = full_path.relative_to(root_dir).as_posix()
            file_list.append((full_path, rel_path))
    return file_list


def create_or_verify_repo(username: str, token: str, repo_name: str, private: bool = False):
    headers = {
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github.v3+json"
    }

    # Check if repo exists
    check_url = f"https://api.github.com/repos/{username}/{repo_name}"
    res = requests.get(check_url, headers=headers)
    if res.status_code == 200:
        print(f"[+] Repository '{username}/{repo_name}' already exists.")
        return res.json()

    # Create repo
    print(f"[*] Creating new repository '{username}/{repo_name}' on GitHub...")
    create_url = "https://api.github.com/user/repos"
    payload = {
        "name": repo_name,
        "description": "Bio-Impact Analyzer - Protein Structural Analysis, AutoDock Vina Docking, & Molecular Dynamics Pipeline",
        "private": private,
        "auto_init": False
    }
    create_res = requests.post(create_url, headers=headers, json=payload)
    if create_res.status_code not in (200, 201):
        raise RuntimeError(f"Failed to create repo: {create_res.status_code} {create_res.text}")

    print(f"[+] Successfully created repository: https://github.com/{username}/{repo_name}")
    return create_res.json()


def upload_repository(username: str, token: str, repo_name: str, commit_message: str = "Initial commit with full Molecular Dynamics & Docking pipeline"):
    headers = {
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github.v3+json"
    }
    base_url = f"https://api.github.com/repos/{username}/{repo_name}"

    # 1. Create or verify repo
    create_or_verify_repo(username, token, repo_name)

    # 2. Collect files
    files = get_files_to_upload(PROJECT_ROOT)
    print(f"[*] Preparing {len(files)} files for upload...")

    # 3. Create blobs for each file
    tree_items = []
    for idx, (full_p, rel_p) in enumerate(files, start=1):
        try:
            with open(full_p, "rb") as f:
                content = f.read()

            # GitHub blob creation
            blob_url = f"{base_url}/git/blobs"
            blob_data = {
                "content": base64.b64encode(content).decode("ascii"),
                "encoding": "base64"
            }
            b_res = requests.post(blob_url, headers=headers, json=blob_data)
            if b_res.status_code in (200, 201):
                sha = b_res.json()["sha"]
                mode = "100755" if full_p.suffix in (".sh", ".exe") or "vina" in rel_p else "100644"
                tree_items.append({
                    "path": rel_p,
                    "mode": mode,
                    "type": "blob",
                    "sha": sha
                })
                print(f"  [{idx}/{len(files)}] Uploaded: {rel_p}")
            else:
                print(f"  [!] Failed blob for {rel_p}: {b_res.status_code} {b_res.text}")
        except Exception as ex:
            print(f"  [!] Error processing {rel_p}: {ex}")

    if not tree_items:
        raise RuntimeError("No files were successfully prepared for upload.")

    # 4. Create git tree
    print("[*] Creating Git Tree...")
    tree_url = f"{base_url}/git/trees"
    tree_res = requests.post(tree_url, headers=headers, json={"tree": tree_items})
    if tree_res.status_code not in (200, 201):
        raise RuntimeError(f"Failed to create git tree: {tree_res.status_code} {tree_res.text}")
    tree_sha = tree_res.json()["sha"]

    # 5. Check if branch exists to get parent commit
    parent_commits = []
    ref_url = f"{base_url}/git/refs/heads/main"
    ref_res = requests.get(ref_url, headers=headers)
    if ref_res.status_code == 200:
        parent_commits.append(ref_res.json()["object"]["sha"])

    # 6. Create commit
    print("[*] Creating Commit...")
    commit_url = f"{base_url}/git/commits"
    commit_payload = {
        "message": commit_message,
        "tree": tree_sha,
        "parents": parent_commits
    }
    commit_res = requests.post(commit_url, headers=headers, json=commit_payload)
    if commit_res.status_code not in (200, 201):
        raise RuntimeError(f"Failed to create commit: {commit_res.status_code} {commit_res.text}")
    commit_sha = commit_res.json()["sha"]

    # 7. Update or create ref (branch heads/main)
    print("[*] Updating branch 'main'...")
    if parent_commits:
        update_ref_url = f"{base_url}/git/refs/heads/main"
        r_res = requests.patch(update_ref_url, headers=headers, json={"sha": commit_sha, "force": True})
    else:
        create_ref_url = f"{base_url}/git/refs"
        r_res = requests.post(create_ref_url, headers=headers, json={"ref": "refs/heads/main", "sha": commit_sha})

    if r_res.status_code not in (200, 201):
        raise RuntimeError(f"Failed to update ref: {r_res.status_code} {r_res.text}")

    repo_url = f"https://github.com/{username}/{repo_name}"
    print(f"\n========================================================")
    print(f"🎉 Successfully pushed project to GitHub!")
    print(f"🔗 Repository URL: {repo_url}")
    print(f"========================================================\n")
    return repo_url


def main():
    parser = argparse.ArgumentParser(description="Upload Bio-Impact project to GitHub.")
    parser.add_argument("--token", help="GitHub Personal Access Token (classic or fine-grained with repo access)")
    parser.add_argument("--username", default="ha9160034-create", help="GitHub username")
    parser.add_argument("--repo", default="Bio-Impact-Analyzer", help="GitHub repository name")
    parser.add_argument("--private", action="store_true", help="Make repository private (default is public)")

    args = parser.parse_args()

    token = args.token or os.environ.get("GITHUB_TOKEN")
    if not token:
        print("\n==========================================================")
        print(" GitHub Authentication Token Required")
        print("==========================================================")
        print("Please provide a Personal Access Token (PAT):")
        print("1. Go to: https://github.com/settings/tokens")
        print("2. Click 'Generate new token (classic)'")
        print("3. Check 'repo' scope and click 'Generate token'")
        print("4. Run: .venv\\Scripts\\python tools\\upload_to_github.py --token YOUR_TOKEN\n")
        try:
            token = input("Enter your GitHub Personal Access Token: ").strip()
        except EOFError:
            pass

    if not token:
        print("Error: No GitHub token provided.")
        sys.exit(1)

    upload_repository(args.username, token, args.repo, private=args.private)


if __name__ == "__main__":
    main()

