console.log(`${Date()}\n${location.origin+location.pathname}index.js loaded`);

(async function(){
    VCF = await import('./export.js')
    //VCF.UI()
    VCF.LDextraction(chrpos1='chr7:16876630',chrpos2='chr7:16863828')
})()