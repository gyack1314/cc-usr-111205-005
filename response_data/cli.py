"""显示原始资料摘要，或执行离散栅格洼地分析。"""
import argparse
import json
from .io import load
from .domain import summarize
from .depression import analyze, to_jsonable


def main():
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=['inspect', 'depressions'])
    p.add_argument('file')
    p.add_argument('--with-paths', action='store_true',
                   help='为每个可排格附无环见证路径（默认仅给等价父结构）')
    a = p.parse_args()
    v, digest = load(a.file)
    if a.command == 'inspect':
        payload = {'source_sha256': digest, 'summary': summarize(v)}
    else:
        result = analyze(v)
        payload = {
            'source_sha256': digest,
            'summary': summarize(v),
            'depression_analysis': to_jsonable(result, with_paths=a.with_paths),
        }
    print(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
