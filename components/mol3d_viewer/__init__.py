"""Custom Streamlit component for 3D molecular visualization with hover screenshot."""

import os
from typing import Optional
import streamlit as st
import streamlit.components.v1 as components

# Declare component pointing to the local directory containing index.html
_component_dir = os.path.dirname(os.path.abspath(__file__))
_component_func = components.declare_component("mol3d_viewer", path=_component_dir)


def render_screenshot_gallery(
    screenshots_key: str = "analysis_screenshots",
    key_prefix: str = "gallery",
) -> None:
    """Render compact thumbnail strip of captured 3D screenshots with delete and clear options.

    Args:
        screenshots_key: Session state key storing the list of captured base64 PNG images.
        key_prefix: Unique key prefix for button instances.
    """
    if screenshots_key not in st.session_state or not isinstance(
        st.session_state[screenshots_key], list
    ):
        st.session_state[screenshots_key] = []

    screenshots = st.session_state[screenshots_key]
    if not screenshots:
        return

    col_hdr_left, col_hdr_right = st.columns([5, 1])
    with col_hdr_left:
        st.caption(f"📸 Captured Views ({len(screenshots)})")
    with col_hdr_right:
        if st.button("Clear All", key=f"{key_prefix}_clear_all", help="Clear all captured screenshots"):
            st.session_state[screenshots_key] = []
            st.rerun()

    # Render thumbnails in rows of up to 6 columns
    max_cols = 6
    to_delete: Optional[int] = None
    for row_start in range(0, len(screenshots), max_cols):
        batch = screenshots[row_start : row_start + max_cols]
        cols = st.columns(max_cols)
        for offset, img in enumerate(batch):
            idx = row_start + offset
            with cols[offset]:
                st.image(img, caption=f"Angle {idx + 1}", use_container_width=True)
                if st.button("🗑️", key=f"{key_prefix}_del_{idx}", help=f"Delete Angle {idx + 1}"):
                    to_delete = idx

    if to_delete is not None:
        st.session_state[screenshots_key].pop(to_delete)
        st.rerun()


def render_3d_viewer(
    viewer_html: str,
    height: int = 500,
    key: str = "mol3d",
    screenshots_key: str = "analysis_screenshots",
    show_gallery: bool = True,
) -> None:
    """Render interactive 3D molecular viewer with hover screenshot capture and optional thumbnail gallery.

    Args:
        viewer_html: HTML content containing the 3Dmol canvas and initialization script.
        height: Height in pixels for the viewer iframe.
        key: Unique key for the Streamlit component instance.
        screenshots_key: Session state key storing the list of captured base64 PNG images.
        show_gallery: If True, renders the screenshot gallery beneath the viewer.
    """
    if not viewer_html:
        return

    # Ensure screenshots list is initialized in session state
    if screenshots_key not in st.session_state or not isinstance(
        st.session_state[screenshots_key], list
    ):
        st.session_state[screenshots_key] = []

    last_ts_key = f"{key}_last_ts"
    if last_ts_key not in st.session_state:
        st.session_state[last_ts_key] = None

    # Render custom component with headless/AppTest fallback
    try:
        component_value = _component_func(
            viewer_html=viewer_html,
            height=height,
            key=key,
            default=None,
        )
    except Exception:
        # In non-standard or headless environments where custom components fail,
        # fail gracefully without crashing the app.
        component_value = None

    # State deduplication: capture screenshot only once per timestamp
    if isinstance(component_value, dict):
        image_data = component_value.get("image")
        ts = component_value.get("ts")
        if image_data and ts and ts != st.session_state.get(last_ts_key):
            st.session_state[last_ts_key] = ts
            st.session_state[screenshots_key].append(image_data)

    if show_gallery:
        render_screenshot_gallery(screenshots_key=screenshots_key, key_prefix=key)


__all__ = ["render_3d_viewer", "render_screenshot_gallery"]
