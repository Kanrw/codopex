"""Single source of truth for the codopex version.

Kept dependency-free so the build backend can read it statically
(``[tool.setuptools.dynamic] version = {attr = "codopex._version.__version__"}``).
"""

__version__ = "0.3.0"
