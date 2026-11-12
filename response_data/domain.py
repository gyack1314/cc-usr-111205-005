"""读取有理数网格，区分屏障与有效高程。"""
from .io import fraction

def summarize(v):
    z=v['elevation'];area=v['area']
    if not z or not z[0] or any(len(r)!=len(z[0]) for r in z) or len(area)!=len(z):raise ValueError('网格形状无效')
    total=0;valid=0
    for row,a in zip(z,area):
        if len(a)!=len(row):raise ValueError('面积形状不符')
        for x,y in zip(row,a):
            if fraction(y)<=0:raise ValueError('格面积应正')
            if x is not None:fraction(x);total+=fraction(y);valid+=1
    for i,j in v['outlets']:
        if not 0<=i<len(z) or not 0<=j<len(z[0]) or z[i][j] is None:raise ValueError('出口不是有效格')
    return {'height':len(z),'width':len(z[0]),'valid_cells':valid,'valid_area':str(total),'outlets':len(v['outlets'])}
