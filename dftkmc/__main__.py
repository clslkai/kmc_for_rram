"""Allow ``python -m dftkmc`` without installation."""
from .cli import main

raise SystemExit(main())
