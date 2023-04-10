console.log('linkage disequilibrium sdk loaded')

var ld = {'vcflib': null}

/**
 * Main global portable module.
 * @namespace 
 * @property {Function} Ld - {@link Ld}
 *
 * @namespace ld
 * @property {Function} calculate_probabilities_snp - {@link ld.calculate_probabilities_snp}
 * @property {Function} calculate_ld - {@link ld.calculate_ld}
 * @property {Function} perform_ld - {@link ld.perform_ld}
 */
 
/** 
* Initializes main library object.
*
* @returns {Object} LD library object.
* 
* @example
* let v = await Ld()
*/
var Ld = async () => {
    var mods = await Promise.all( [import('https://episphere.github.io/vcf/export.js'), loadScript('https://cdn.jsdelivr.net/npm/jstat@latest/dist/jstat.min.js')] )
    ld.vcflib = mods[0]
}

/** 
* Calculates probabilities and counts occurrences of each phenotype
*
* @param {Array} snp_info Array containing a line of the vcf file concerning the snp
*
* @returns {Object} Object containing information of count (homozygous dominant, heterozygous and homozygous recessive) and probabilities p and q for a certain snp.
* 
* @example
* let v = await ld.calculate_probabilities_snp()
*/
ld.calculate_probabilities_snp = (snp_info) => {
    var probs={}
    var j = 0
    var index_gen = 0
    snp_info.forEach( e => { if(e.length==3 && e.indexOf('|')==1 && index_gen==0){ index_gen=j; } j+=1  })
    if(index_gen!=0){
        var cnts={}
        var total = 0
        j=0
        snp_info.forEach( e => { 
            if(e.length==3 && e.indexOf('|')==1 && j>=index_gen ){ 
                if( e.split('|')[0] != e.split('|')[1] ){
                    e='mix'
                }
                if( ! Object.keys(cnts).includes(e) ){
                    cnts[e]=0
                }
                cnts[e]+=1 
                total+=1
            } 
            j+=1  
        })
        total *= 2 // Total number of chromosomes in the population
        
        for ( var k of ['0|0', 'mix', '1|1']){
            if( ! Object.keys(cnts).includes(k) ){
                cnts[k]=0
            }
        }
        
        for(var k of Object.keys(cnts)){
            probs['count_'+k] = cnts[k]
        }
        probs['p'] = ( (2*cnts['0|0']) + cnts['mix'] )/total
        probs['q'] = ( cnts['mix'] + (2*cnts['1|1']) )/total
    }
    return probs
}

/** 
* Calculates r2 and d' metrics for a pair of snps
*
* @param {Array} snp1 Array containing a line of the vcf file concerning the first snp
* @param {Array} snp2 Array containing a line of the vcf file concerning the second snp
*
* @returns {Object} Object containing information of haplotypes count and the following statistical measures: chisquare, D', r², D and p-value.
* 
* @example
* let v = calculald.te_ld()
*/
ld.calculate_ld = (snp1, snp2) => {
    var result = {'dl': -1, 'r2': -1}
    
    var index_gen = 0
    var j=0
    snp1.forEach( e => { if(e.length==3 && e.indexOf('|')==1 && index_gen==0){ index_gen=j; } j+=1  })
    if(index_gen!=0 && snp1.length == snp2.length){
        var cnts={}
        j=0
        var snp = snp1
        var snp_2=snp2
        if(snp1.length < snp2.length){
            snp=snp2
            snp_2=snp1
        }
        snp.forEach( e => { 
            if(e.length==3 && e.indexOf('|')==1 && j>=index_gen ){ 
                if( e.split('|')[0] != e.split('|')[1] ){
                    e='mix'
                }
                var f = snp_2[j]
                if(f!=undefined){
                    if( f.split('|')[0] != f.split('|')[1] ){
                        f='mix'
                    }
                    
                    var aux = []
                    
                    if(e=='0|0'){
                        if(f=='mix'){
                            aux.push('00')
                            aux.push('01')
                        }
                        if(f=='0|0'){
                            aux.push('00')
                        }
                        if(f=='1|1'){
                            aux.push('01')
                        }
                    } 
                    if(e=='1|1' ){
                        if(f=='mix'){
                            aux.push('10')
                            aux.push('11')
                        }
                        if(f=='0|0'){
                            aux.push('10')
                        }
                        if(f=='1|1'){
                            aux.push('11')
                        }
                    } 
                    if(e=='mix' ){
                        if(f=='mix'){
                            aux.push('00')
                            aux.push('01')
                            aux.push('10')
                            aux.push('11')
                        }
                        if(f=='0|0'){
                            aux.push('00')
                            aux.push('10')
                        }
                        if(f=='1|1'){
                            aux.push('01')
                            aux.push('11')
                        }
                    }
                    
                    for (var x of aux){
                        if( ! Object.keys(cnts).includes(x) ){
                            cnts[x]=0
                        }
                        cnts[x]+=1 
                    }
                }
            } 
            
            j+=1  
        })
        j*=2
        
        var probs1 = ld.calculate_probabilities_snp(snp1)
        var probs2 = ld.calculate_probabilities_snp(snp2)
        
        var result = {}
        result['chisq'] = 0
        for( var k of Object.keys(cnts) ){
            result['p'+k]=cnts[k]/j
            
            var p1=0
            var p2=0
            if(k[0]=='0'){
                p1=probs1['p']
            }
            if(k[0]=='1'){
                p1=probs1['q']
            }
            if(k[1]=='0'){
                p2=probs2['p']
            }
            if(k[1]=='1'){
                p2=probs2['q']
            }
            var expected = p1*p2*j
            console.log(k, expected)
            result['chisq'] += Math.pow( (cnts[k]-expected), 2)/expected
        }
        result['pvalue'] = jStat.chisquare.pdf( result['chisq'], 1 )
        
        result['d'] = (result['p00']*result['p11']) - (result['p01']*result['p10'])
        var denominator = Math.min(probs1['p']*probs2['q'], probs1['q']*probs2['p'])  // dmax
        if( result['d']<0 ){
            denominator = Math.max(probs1['p']*probs2['q'], probs1['q']*probs2['p']) // dmin
        }
        result['dl'] = result['d']/denominator
        var denominator_r2 = probs1['p']*probs1['q']*probs2['p']*probs2['q']
        result['r2'] = result['d']/denominator_r2
    }
    
    return result
}

/** 
* Perform pairwise LD analysis with the snps of a certain chromosome in a given position range
*
* @param {number} chromosome Chromosome
* @param {number} start Start Position
* @param {number} end End Position
*
* @returns {Object} Object containing the linkage disequilibium result for all pairs combination of SNPs found in the given posiion range.
* 
* @example
* let v = await ld.perform_ld()
*/
ld.perform_ld = (chrom, start, end) => {
    var chrom = ld_chrom.value
    var start = ld_start.value
    var end = ld_end.value
    
    var ld_result={}

    var url=`http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
    if(chrom=='MT'){
        url='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz'
    }
    if(chrom=='Y'){
        url='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz'
    }
    if(chrom=='X'){
        url='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz'
    }
    
    if(start < end){
        action_ld.disabled=true
        action_ld.innerHTML='Analyzing ...'
        console.log(ld.vcflib)
        ld.vcflib.Vcf(url).then( async (value) => {
            var vld = value
            var queries=[]
            for( var i=start; i<=end; i++){
                queries.push(chrom+','+i)
            }
            var snps={}
            var info = await Promise.all( queries.map( async e => {
                var r = await vld.query(e)
                //console.log(e, '----', r)
                if(r!=undefined){
                    if(r.hit.length>0){
                        r=r.hit[0]
                        var ide = r[2]
                        if(r[2]=='.'){
                            ide=r[0]+'_'+r[1]
                        }
                        snps[ide] = r
                        return ide
                    }
                }
                
                return null
            } ) )
            console.log(snps)
            console.log(info)
            var keys = Object.keys(snps)
            var i=0
            for (var k of keys){
                var j=0
                for (var v of keys){
                    if(i<j){
                        ld_result[k+'-'+v] = ld.calculate_ld(snps[k], snps[v])
                    }
                    j+=1
                }
                i+=1
            }
            
            action_ld.disabled=false
            action_ld.innerHTML='Analyze'
            
            console.log(ld_result)
            
            return ld_result
        })
    }
    else{
        alert('End position must be higher than the start position')
        return ld_result
    }
    
}

/** 
* Load a certain dependency library from link
* 
*
* @param {string} url Library URL.
* 
* @example
* loadScript('https://cdnjs.cloudflare.com/ajax/libs/pako/1.0.11/pako.min.js')
*
*/
ld.loadScript= async function(url){
	console.log(`${url} loaded`)
    async function asyncScript(url){
        let load = new Promise((resolve,regect)=>{
            let s = document.createElement('script')
            s.src=url
            s.onload=resolve
            document.head.appendChild(s)
        })
        await load
    }
    // satisfy dependencies
    await asyncScript(url)
} 

if(typeof(jStat)=="undefined"){
    ld.loadScript('https://cdn.jsdelivr.net/npm/jstat@latest/dist/jstat.min.js')
}

if(ld.vcflib==null){
    import('https://episphere.github.io/vcf/export.js').then( mod => {
        ld.vcflib = mod
    })
}

