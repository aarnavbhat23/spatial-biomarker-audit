"""A local HTTP observation service, with simulated costs and hidden reference calls."""
import json,threading,time
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.request import Request,urlopen
import numpy as np
import simpy
from joint_review import choose_action,intervals


class ObservationService:
    def __init__(self,weights,totals,owners,recorded,reference,budget,review_cost=.05,feedback=.1):
        self.weights=weights;self.totals=totals;self.owners=owners;self._recorded=recorded;self._reference=reference
        self.budget=budget;self.review_cost=review_cost;self.feedback=feedback
        self.seen=np.zeros(len(owners),bool);self.reviewed=self.seen.copy();self.env=simpy.Environment()
        self.scanner=simpy.Resource(self.env,capacity=1);self.reviewer=simpy.Resource(self.env,capacity=1)
        self.events=[];self._lock=threading.Lock()
    def metadata(self):
        return {'weights':self.weights.tolist(),'totals':self.totals.tolist(),'owners':self.owners.tolist(),
                'budget':self.budget,'review_cost':self.review_cost,'field_cost':1+self.feedback,'cost_units':'hypothetical field equivalents',
                'baseline':'B calls and graph already known exactly; their common acquisition cost is excluded'}
    def observe(self,kind,index):
        with self._lock:return self._observe(kind,index)
    def _observe(self,kind,index):
        start=self.env.now
        if kind=='field':
            ids=np.flatnonzero(self.owners==index)
            new=ids[~self.seen[ids]]
            if not np.any(self.owners==index):raise ValueError('Unknown field')
            duration=(1+self.feedback) if len(new) else 0.;resource=self.scanner
            self.seen[ids]=True;values=self._recorded[ids]
        elif kind=='review':
            if index<0 or index>=len(self.seen) or not self.seen[index]:raise ValueError('Review requires an acquired cell')
            ids=np.array([index]);duration=0. if self.reviewed[index] else self.review_cost;resource=self.reviewer
            self.reviewed[index]=True;values=self._reference[ids]
        else:raise ValueError('Unknown observation type')
        def job():
            with resource.request() as request:
                yield request;yield self.env.timeout(duration)
        process=self.env.process(job());self.env.run(until=process)
        event={'kind':kind,'index':int(index),'ids':ids.tolist(),'calls':values.tolist(),
               'simulated_start':start,'simulated_end':self.env.now,'simulated_increment':duration}
        self.events.append(event);return event


@contextmanager
def http_service(service):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def reply(self,status,data):
            raw=json.dumps(data).encode();self.send_response(status);self.send_header('Content-Type','application/json')
            self.send_header('Content-Length',str(len(raw)));self.end_headers();self.wfile.write(raw)
        def do_GET(self):
            self.reply(200,service.metadata()) if self.path=='/metadata' else self.reply(404,{'error':'Unknown route'})
        def do_POST(self):
            try:
                body=json.loads(self.rfile.read(int(self.headers['Content-Length'])))
                if self.path not in ('/field','/review'):raise ValueError('Unknown route')
                result=service.observe(self.path[1:],int(body['id']));self.reply(200,result)
            except (ValueError,KeyError,IndexError) as exc:self.reply(400,{'error':str(exc)})
    server=ThreadingHTTPServer(('127.0.0.1',0),Handler);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    try:yield 'http://127.0.0.1:'+str(server.server_port)
    finally:server.shutdown();server.server_close();thread.join()


class HttpClient:
    def __init__(self,url):self.url=url
    def metadata(self):
        with urlopen(self.url+'/metadata',timeout=10) as r:return json.load(r)
    def observe(self,kind,index):
        request=Request(self.url+'/'+kind,data=json.dumps({'id':int(index)}).encode(),headers={'Content-Type':'application/json'})
        with urlopen(request,timeout=10) as r:return json.load(r)


def replay(client,policy,width=.005,seed=0,max_actions=10000):
    meta=client.metadata();weights=np.array(meta['weights']);totals=np.array(meta['totals']);owners=np.array(meta['owners'])
    n=len(owners);seen=np.zeros(n,bool);reviewed=seen.copy();recorded=seen.copy();verified=seen.copy()
    rng=np.random.default_rng(seed);checkpoints=[];start=time.perf_counter()
    current=intervals(weights,totals,recorded,seen,reviewed,verified,meta['budget'])
    for _ in range(max_actions):
        if np.max(np.diff(current,axis=1))<=width:break
        action=choose_action(weights,totals,owners,recorded,seen,reviewed,verified,meta['budget'],meta['review_cost'],policy,rng,meta['field_cost'])
        if action is None:break
        kind,index=action;event=client.observe(kind,index);ids=event['ids']
        if kind=='field':seen[ids]=True;recorded[ids]=event['calls']
        else:reviewed[ids]=True;verified[ids]=event['calls']
        current=intervals(weights,totals,recorded,seen,reviewed,verified,meta['budget'])
        checkpoints.append({'kind':kind,'index':index,'intervals':current.tolist(),'simulated_time':event['simulated_end']})
    return {'policy':policy,'fields':len({c['index'] for c in checkpoints if c['kind']=='field'}),'reviews':int(reviewed.sum()),
            'simulated_time':checkpoints[-1]['simulated_time'] if checkpoints else 0.,'software_wall_seconds':time.perf_counter()-start,
            'intervals':current.tolist(),'stopped':bool(np.max(np.diff(current,axis=1))<=width),'checkpoints':checkpoints}
