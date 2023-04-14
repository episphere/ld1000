console.log(`${Date()}\n${location.origin+location.pathname}index.js loaded`);

(async function(){
    V = await import('./export.js')
    V.UI()
})()