"""Unit tests for mol3d_viewer custom component."""

import unittest
from unittest.mock import patch
import streamlit as st
from streamlit.testing.v1 import AppTest

import components.mol3d_viewer as mol3d_viewer


class TestMol3DViewer(unittest.TestCase):
    """Test suite for mol3d_viewer Streamlit component."""

    def test_empty_viewer_html(self):
        """Verify render_3d_viewer returns without error on empty HTML."""
        mol3d_viewer.render_3d_viewer("", key="empty_test")
        mol3d_viewer.render_3d_viewer(None, key="none_test")

    def test_apptest_headless_mode(self):
        """Verify render_3d_viewer runs cleanly in headless AppTest mode."""
        script = """
import streamlit as st
from components.mol3d_viewer import render_3d_viewer
render_3d_viewer('<div>test_mol</div>', key='apptest_key', screenshots_key='apptest_shots')
"""
        at = AppTest.from_string(script).run()
        self.assertEqual(len(at.exception), 0)
        self.assertIn("apptest_shots", at.session_state)
        self.assertEqual(at.session_state["apptest_shots"], [])

    def test_thumbnail_strip_and_captions(self):
        """Verify thumbnail strip renders images, captions, and buttons."""
        script = """
import streamlit as st
from components.mol3d_viewer import render_3d_viewer
if 'shots' not in st.session_state:
    st.session_state['shots'] = ['data:image/png;base64,aaa', 'data:image/png;base64,bbb']
render_3d_viewer('<div>test_mol</div>', key='strip_key', screenshots_key='shots')
"""
        at = AppTest.from_string(script).run()
        self.assertEqual(len(at.exception), 0)
        self.assertEqual(len(at.image), 2)
        self.assertEqual(at.image[0].captions, ["Angle 1"])
        self.assertEqual(at.image[1].captions, ["Angle 2"])
        # Buttons: 1 Clear All + 2 individual deletes
        self.assertEqual(len(at.button), 3)

    def test_clear_all_button(self):
        """Verify Clear All button empties screenshots list."""
        script = """
import streamlit as st
from components.mol3d_viewer import render_3d_viewer
if 'shots' not in st.session_state:
    st.session_state['shots'] = ['data:image/png;base64,aaa', 'data:image/png;base64,bbb']
render_3d_viewer('<div>test_mol</div>', key='strip_key', screenshots_key='shots')
"""
        at = AppTest.from_string(script).run()
        # Click Clear All (first button)
        at.button[0].click().run()
        self.assertEqual(len(at.exception), 0)
        self.assertEqual(at.session_state["shots"], [])
        self.assertEqual(len(at.image), 0)

    def test_delete_individual_thumbnail(self):
        """Verify individual delete button removes only that screenshot."""
        script = """
import streamlit as st
from components.mol3d_viewer import render_3d_viewer
if 'shots' not in st.session_state:
    st.session_state['shots'] = ['data:image/png;base64,111', 'data:image/png;base64,222']
render_3d_viewer('<div>test_mol</div>', key='strip_key', screenshots_key='shots')
"""
        at = AppTest.from_string(script).run()
        # button 0 is Clear All, button 1 is delete first item
        at.button[1].click().run()
        self.assertEqual(len(at.exception), 0)
        self.assertEqual(at.session_state["shots"], ["data:image/png;base64,222"])

    def test_state_deduplication(self):
        """Verify timestamp deduplication avoids duplicating captured screenshots."""
        st.session_state.clear()
        with patch("components.mol3d_viewer._component_func") as mock_func:
            # 1. First capture event
            mock_func.return_value = {"image": "data:image/png;base64,shotA", "ts": 1000}
            mol3d_viewer.render_3d_viewer("<div>mol</div>", key="k_dedup", screenshots_key="s_dedup")
            self.assertEqual(st.session_state["s_dedup"], ["data:image/png;base64,shotA"])
            self.assertEqual(st.session_state["k_dedup_last_ts"], 1000)

            # 2. Rerun with identical ts does not duplicate
            mol3d_viewer.render_3d_viewer("<div>mol</div>", key="k_dedup", screenshots_key="s_dedup")
            self.assertEqual(len(st.session_state["s_dedup"]), 1)

            # 3. Subsequent capture event with newer ts appends
            mock_func.return_value = {"image": "data:image/png;base64,shotB", "ts": 2000}
            mol3d_viewer.render_3d_viewer("<div>mol</div>", key="k_dedup", screenshots_key="s_dedup")
            self.assertEqual(len(st.session_state["s_dedup"]), 2)
            self.assertEqual(st.session_state["s_dedup"][1], "data:image/png;base64,shotB")
            self.assertEqual(st.session_state["k_dedup_last_ts"], 2000)

    def test_graceful_component_failure(self):
        """Verify exceptions from _component_func are swallowed gracefully."""
        st.session_state.clear()
        with patch("components.mol3d_viewer._component_func", side_effect=RuntimeError("Component crashed")):
            # Should not raise
            mol3d_viewer.render_3d_viewer("<div>mol</div>", key="k_fail", screenshots_key="s_fail")
            self.assertEqual(st.session_state["s_fail"], [])

    def test_show_gallery_false_and_standalone_gallery(self):
        """Verify show_gallery=False suppresses in-viewer gallery, and render_screenshot_gallery renders standalone."""
        script = """
import streamlit as st
from components.mol3d_viewer import render_3d_viewer, render_screenshot_gallery
if 'shots' not in st.session_state:
    st.session_state['shots'] = ['data:image/png;base64,aaa', 'data:image/png;base64,bbb']
# Viewer without gallery
render_3d_viewer('<div>test_mol</div>', key='view_no_gal', screenshots_key='shots', show_gallery=False)
# Separate unified gallery
render_screenshot_gallery(screenshots_key='shots', key_prefix='unified')
"""
        at = AppTest.from_string(script).run()
        self.assertEqual(len(at.exception), 0)
        # Images rendered by standalone gallery
        self.assertEqual(len(at.image), 2)
        self.assertEqual(at.image[0].captions, ["Angle 1"])
        self.assertEqual(at.image[1].captions, ["Angle 2"])
        # Buttons: 1 Clear All + 2 deletes
        self.assertEqual(len(at.button), 3)
        # Clear All button key has prefix 'unified'
        self.assertIn("unified_clear_all", at.button[0].key)


if __name__ == "__main__":
    unittest.main()
