"""Execute invented-data unit tests only."""
from pathlib import Path
import os,subprocess,sys
if __name__=='__main__':
 root=Path(__file__).resolve().parent
 raise SystemExit(subprocess.run([sys.executable,'-B',str(root/'test_paired_stats_v1.py'),'-v'],
    cwd=root,timeout=120,check=False,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'}).returncode)
