import {loadPyodide} from 'pyodide';
import fs from 'node:fs';
const p=await loadPyodide();
const base=new URL('../src/dynacomp/',import.meta.url).pathname;p.FS.mkdirTree('/code/dynacomp');
for(const name of fs.readdirSync(base))if(name.endsWith('.py'))p.FS.writeFile('/code/dynacomp/'+name,fs.readFileSync(base+'/'+name,'utf8'));
const result=await p.runPythonAsync(`import sys,types,json\npkg=types.ModuleType('dynacomp');pkg.__path__=['/code/dynacomp'];sys.modules['dynacomp']=pkg\nfrom dynacomp.explorer import scenario\nresults=[scenario({'example':kind,'output':out,'steps':8}) for kind in ['warmth','supplies'] for out in ['whole','one']]\nassert all(r['full']==r['reduced'] for r in results)\njson.dumps([{'title':r['title'],'retained':r['retained_values'],'max_difference':r['max_actual_difference']} for r in results])`);
console.log(result);
