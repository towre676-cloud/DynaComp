import {loadPyodide} from 'https://cdn.jsdelivr.net/pyodide/v0.29.5/full/pyodide.mjs';
let runtime;
async function prepare(bundle){
 if(runtime)return runtime;
 self.postMessage({status:'preparing'});
 const py=await loadPyodide();py.FS.mkdirTree('/model/dynacomp');
 for(const [name,source] of Object.entries(bundle.engine_files)){
  if(!/^dynacomp\/[a-z_]+\.py$/.test(name))throw Error('Unsupported engine file path');
  const digest=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(source))),x=>x.toString(16).padStart(2,'0')).join('');
  if(digest!==bundle.source_hashes[name])throw Error('Engine source did not match its public digest');
  py.FS.writeFile('/model/'+name,source);
 }
 // Activate only the portable package namespace. Desktop native-library
 // initialization is not executed in WebAssembly; the mathematical modules
 // below remain the original source-bound DynaComp implementations.
 await py.runPythonAsync("import sys,types,json\npkg=types.ModuleType('dynacomp');pkg.__path__=['/model/dynacomp'];sys.modules['dynacomp']=pkg\nfrom dynacomp.explorer import scenario");
 runtime=py;return py;
}
self.onmessage=async event=>{try{const py=await prepare(event.data.bundle);py.globals.set('request_json',JSON.stringify(event.data.request));const raw=await py.runPythonAsync('json.dumps(scenario(json.loads(request_json)))');self.postMessage({status:'complete',result:JSON.parse(raw)});}catch(e){self.postMessage({status:'error',message:String(e)});}};
