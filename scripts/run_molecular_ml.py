"""Launch local CPU experiments; expose the installed PyTorch OpenMP library on macOS."""
import importlib.util,os,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
env=dict(os.environ,PYTHONPATH=str(root/'src'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
if sys.platform=='darwin':
    spec=importlib.util.find_spec('torch')
    lib=Path(spec.origin).parent/'lib'
    env['DYLD_LIBRARY_PATH']=str(lib)+(os.pathsep+env['DYLD_LIBRARY_PATH'] if env.get('DYLD_LIBRARY_PATH') else '')
os.execve(sys.executable,[sys.executable,'-m','pktk.ml.run',*sys.argv[1:]],env)
