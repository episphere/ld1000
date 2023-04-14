import {Vcf} from 'https://episphere.github.io/vcf/export.js'

async function connectVCF(url=1){
    const baseURL='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
    if(typeof(url)=='number'){
        url=baseURL.replace('chr1',`chr${url}`)
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
        Chromossome 1: <select id="chrSelect1"></select> Position 1 <input type="number" id="pos1" value=10177>
        <br>Chromossome 2<sup>*</sup>: <select id="chrSelect2"></select>  Position 2 <input type="text" id="pos2" value=10352>
        <p><button>Retrieve positions</button> directly from <a href="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" target="_blank">1000genomes</a>:</p>
        <div>...</div>
        <p>* <span style="font-size:x-small">by default the same chromossome, different chromossomes are not linked.</span></p>`
    let chrs=[...Array(22)].map((_,i)=>`chr${i+1}`).concat(['X','Y','MT'])
    let chrSelect1=div.querySelector('#chrSelect1');
    chrs.forEach(chr=>{
        let opt = document.createElement('option')
        opt.textContent=chr
        chrSelect1.appendChild(opt)
    })
    let chrSelect2=div.querySelector('#chrSelect2');
    chrs.forEach(chr=>{
        let opt = document.createElement('option')
        opt.textContent=chr
        chrSelect2.appendChild(opt)
    })
    chrSelect1.onchange=function(){
         chrSelect2.value=chrSelect1.value
    }
}

export{
    connectVCF,
    UI
}