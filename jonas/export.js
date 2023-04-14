import {Vcf} from 'https://episphere.github.io/vcf/export.js'

var VV={} // Caching connections here

async function connectVCF(url=1){
    const baseURL='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
    if(typeof(url)=='number'){
        url=baseURL.replace('chr1',`chr${url}`)
    }else if((typeof(url)=='string')&(url.length<6)){
        url=baseURL.replace('chr1',`${url}`)
    }
    if(url.length<3){  // MT X Y
        switch (url){
            case 'MT':
                url='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz';
                break;
            case 'X':
                url='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz';
                break;
            case 'Y':
                url='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz';
                break;
            default:
                error(`${url} not found`)
        }
    }
    return await Vcf(url)
}

async function UI(div){
    if(!div){
        div = document.createElement('div')
        div.id = 'LD1000UI'
        document.body.appendChild(div)
    }
    if(typeof(div)=='string'){
        div = document.getElementById(div)
    }
    
    div.innerHTML=`
        <h3><a href="https://episphere.github.io/ld1000/jonas" target="_blank">LD1000</a> calculator</h3>
        Chromossome 1: <select id="chrSelect1"></select> Position 1 <input type="number" id="pos1" value=16876630>
        <br>Chromossome 2<sup>*</sup>: <select id="chrSelect2"></select>  Position 2 <input type="text" id="pos2" value=16863828>
        <p><button id="retrievePosButton">Retrieve positions</button> directly from <a href="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" target="_blank">1000genomes</a>:</p>
        <div id="resultsDiv">...</div>
        <p>* <span style="font-size:x-small">by default the same chromossome, different chromossomes are not linked.</span></p>`
    let chrs=[...Array(22)].map((_,i)=>`chr${i+1}`).concat(['X','Y','MT'])
    let chrSelect1=div.querySelector('#chrSelect1');
    chrs.forEach(chr=>{
        let opt = document.createElement('option')
        opt.textContent=chr
        chrSelect1.appendChild(opt)
    })
    chrSelect1.value='chr7'
    let chrSelect2=div.querySelector('#chrSelect2');
    chrs.forEach(chr=>{
        let opt = document.createElement('option')
        opt.textContent=chr
        chrSelect2.appendChild(opt)
    })
    chrSelect2.value='chr7'
    chrSelect1.onchange=function(){
         chrSelect2.value=chrSelect1.value
    }

    async function retrieveChrPos(){ // caching VCF connections
        if(!VV[chrSelect1.value]){
            VV[chrSelect1.value] = await connectVCF(chrSelect1.value)
        }
        if(!VV[chrSelect2.value]){ // in the rare instance chrSelect2 is different from 1
            VV[chrSelect2.value] = await connectVCF(chrSelect2.value)
        }
        let V1 = VV[chrSelect1.value]
        let V2 = VV[chrSelect2.value]
        let pos1 = div.querySelector('#pos1').value
        let pos2 = div.querySelector('#pos2').value
        let chr1 = chrSelect1.value.match(/[0-9,X,Y,M,T]+/)[0]
        let chr2 = chrSelect2.value.match(/[0-9,X,Y,M,T]+/)[0]
        let q1 = await V1.query(`${chr1}:${pos1}`)
        let q2 = await V2.query(`${chr2}:${pos2}`)
        debugger
        div.querySelector('#resultsDiv').innerHTML=`<hr>${JSON.stringify(q1.hit)}<hr>${JSON.stringify(q2.hit)}<hr>`
        // try 7:16876630 vs 7:16876630

        
    }

    div.querySelector('#retrievePosButton').onclick=retrieveChrPos

    
    //retrieveChrPos()
}

export{
    connectVCF,
    UI,
    VV
}