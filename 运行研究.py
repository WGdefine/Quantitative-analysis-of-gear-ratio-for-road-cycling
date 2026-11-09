"""Find Python 3.10+ and run this study. This launcher also works on Python 3.4.

Examples: python 运行研究.py --mode All
          python 运行研究.py --mode Check
Set BICYCLE_PYTHON to select an interpreter explicitly.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['All', 'Check', 'Test', 'Data', 'Figures'], default='All')
    args = parser.parse_args()
    root = os.path.abspath(os.path.dirname(__file__))
    candidates = []
    override = os.environ.get('BICYCLE_PYTHON')
    if override:
        candidates.append(override)
    candidates.extend([
        os.path.join(os.path.expanduser('~'), '.cache', 'codex-runtimes',
                     'codex-primary-runtime', 'dependencies', 'python', 'python.exe'),
        sys.executable,
    ])
    candidates.extend(glob.glob(os.path.join(os.environ.get('LOCALAPPDATA', ''), 'Programs', 'Python', 'Python*', 'python.exe')))
    candidates.extend([shutil.which('python3'), shutil.which('python')])
    probe = 'import sys; assert sys.version_info >= (3,10)'
    if args.mode in ('All', 'Figures'):
        probe += '; import matplotlib.pyplot'
    selected, seen = None, set()
    for candidate in candidates:
        if not candidate or candidate in seen or not os.path.isfile(candidate):
            continue
        seen.add(candidate)
        try:
            result = subprocess.call([candidate, '-c', probe], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except OSError:
            continue
        if result == 0:
            selected = candidate
            break
    if selected is None:
        print('Python 3.10+ is required; All/Figures also require matplotlib. Set BICYCLE_PYTHON to its executable.')
        return 1
    print('Python: ' + selected)
    sys.stdout.flush()
    tasks = {
        'Check': [['复算模型.py', '--check']],
        'Test': [['-m', 'unittest', '-v']],
        'Data': [['复算模型.py']],
        'Figures': [['绘制图表.py']],
        'All': [['-m', 'unittest', '-v'], ['更新发布稿.py']],
    }
    for task in tasks[args.mode]:
        result = subprocess.call([selected, '-B', '-X', 'utf8'] + task, cwd=root)
        if result:
            return result
    return 0


if __name__ == '__main__':
    sys.exit(main())
