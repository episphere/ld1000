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

export{
    connectVCF
}