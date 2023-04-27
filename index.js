console.log(`${Date()}\n${location.origin+location.pathname}index.js loaded`);

(async function(){
    VCF = await import('./export.js')
    VCF.UI()
    //data = await VCF.LDextraction()
})()