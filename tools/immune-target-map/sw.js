"use strict";

const CACHE_PREFIX = "immune-target-map-";
const CACHE = CACHE_PREFIX + "v1";
const ASSETS = [
  "./index.html",
  "./manifest.json",
  "./icon.svg",
  "./icon-180.png",
  "./icon-192.png",
  "./icon-512.png",
  "./icon-maskable-512.png"
];
const scope = new URL(self.registration.scope);
const inScope = url => url.origin === scope.origin && url.pathname.startsWith(scope.pathname);

self.addEventListener("install", event => {
  event.waitUntil(
    caches.open(CACHE).then(cache => cache.addAll(ASSETS)).then(() => self.skipWaiting())
  );
});

self.addEventListener("activate", event => {
  event.waitUntil(
    caches.keys()
      .then(keys => Promise.all(
        keys
          .filter(key => key.startsWith(CACHE_PREFIX) && key !== CACHE)
          .map(key => caches.delete(key))
      ))
      .then(() => self.clients.claim())
  );
});

self.addEventListener("fetch", event => {
  const request = event.request;
  if (request.method !== "GET") return;

  const url = new URL(request.url);
  if (!inScope(url)) return;

  if (request.mode === "navigate") {
    event.respondWith(
      fetch(request)
        .then(response => {
          if (response.ok) {
            const copy = response.clone();
            event.waitUntil(caches.open(CACHE).then(cache => cache.put("./index.html", copy)));
          }
          return response;
        })
        .catch(() => caches.match("./index.html").then(hit => hit || Response.error()))
    );
    return;
  }

  event.respondWith(
    caches.open(CACHE).then(cache =>
      cache.match(request).then(hit => {
        const network = fetch(request)
          .then(response => {
            if (response.ok && response.type === "basic") {
              return cache.put(request, response.clone()).then(() => response);
            }
            return response;
          })
          .catch(() => null);
        if (hit) {
          event.waitUntil(network);
          return hit;
        }
        return network.then(response => response || Response.error());
      })
    )
  );
});
