# conftest.py --- pytest session setup shared by all tests under tests/.
#
# Force a non-interactive matplotlib backend so that any plt.show() -- whether
# in a test or inside library code such as rdesigneur.display()
# (rdesigneur.py: plt.show(block=...)) -- becomes a no-op instead of popping a
# blocking window and hanging the test run. Set both the env var (honoured
# before matplotlib is first imported) and use() (in case it already was).
import os

os.environ.setdefault("MPLBACKEND", "Agg")

try:
    import matplotlib
    matplotlib.use("Agg", force=True)
except Exception:
    pass
