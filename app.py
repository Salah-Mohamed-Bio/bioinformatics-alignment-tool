"""
Minimal wrapper - imports the Flask app from api/index.py.
This allows Vercel's auto-detection to find 'app' at the top level.
"""
import sys
import os

# Ensure the project root and api directory are on the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from api.index import app as app