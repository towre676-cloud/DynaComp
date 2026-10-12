const CACHE='dynacomp-shell-1';
self.addEventListener('install',e=>e.waitUntil(caches.open(CACHE).then(c=>c.addAll(['/','/style.css','/app.js','/engine-worker.mjs','/icon.svg','/manifest.webmanifest']))));
self.addEventListener('activate',e=>e.waitUntil(self.clients.claim()));
self.addEventListener('fetch',e=>{if(e.request.method!=='GET')return;const u=new URL(e.request.url);if(u.origin!==self.location.origin&&!u.hostname.endsWith('jsdelivr.net'))return;e.respondWith((async()=>{const cache=await caches.open(CACHE);try{const response=await fetch(e.request);if(response.ok)await cache.put(e.request,response.clone());return response;}catch(error){const cached=await cache.match(e.request);if(cached)return cached;throw error;}})());});
