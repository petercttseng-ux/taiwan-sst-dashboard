import os,sys,re,json,time
import numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sst_core as C
ROOT=os.environ.get('SST_ROOT','./衛星水溫圖')
OUT=os.environ.get('SST_OUT','./out')
os.makedirs(OUT,exist_ok=True)
DATE=re.compile(r'(20\d{2})[-_]?(\d{2})[-_]?(\d{2})')
def collect():
    it={}
    for d in sorted(os.listdir(ROOT)):
        p=os.path.join(ROOT,d)
        if not os.path.isdir(p): continue
        for f in sorted(os.listdir(p)):
            fl=f.lower()
            if not fl.endswith('.png') or 'contour_only' in fl or 'onlly' in fl: continue
            if d[0]=='S':
                if 'domain-taiwan' not in f: continue
                m=DATE.search(d)
            else:
                m=DATE.search(f)
            if not m: continue
            k='%s-%s-%s'%m.groups()
            fp=os.path.join(p,f); sz=os.path.getsize(fp)
            if k not in it or sz>it[k][1]: it[k]=(fp,sz)
    return sorted(it.items())
logp=os.path.join(OUT,'log.jsonl')
done=set()
if os.path.exists(logp):
    for l in open(logp,encoding='utf-8'):
        try: done.add(json.loads(l)['date'])
        except: pass
cube={}
cp=os.path.join(OUT,'cube.npz')
if os.path.exists(cp):
    z=np.load(cp)
    for k in z.files: cube[k]=z[k]
todo=[(k,v) for k,v in collect() if k not in done or k not in cube]
print('todo',len(todo),flush=True)
log=open(logp,'a',encoding='utf-8'); t0=time.time(); n=0
for date,(fp,sz) in todo:
    try:
        T,info=C.decode(fp); cube[date]=T.astype(np.float32)
        rec=dict(date=date,ok=True,file=os.path.basename(fp),**info)
    except Exception as e:
        rec=dict(date=date,ok=False,file=os.path.basename(fp),err=str(e)[:150])
    log.write(json.dumps(rec,ensure_ascii=False)+'\n'); log.flush(); n+=1
    if n%50==0:
        el=time.time()-t0
        print('%d/%d %.1f/s'%(n,len(todo),n/el),flush=True)
    if n%300==0: np.savez_compressed(cp,**cube)
np.savez_compressed(cp,**cube)
print('DONE n=%d cube=%d %.0fs'%(n,len(cube),time.time()-t0),flush=True)
