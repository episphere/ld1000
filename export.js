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

async function LDextraction(chrpos1='chr7:16876630',chrpos2='chr7:16863828'){
    let chr1 = chrpos1.match(/(chr\w+):(\w+)/)[1]
    let pos1 = chrpos1.match(/(chr\w+):(\w+)/)[2]
    let chr2 = chrpos2.match(/(chr\w+):(\w+)/)[1]
    let pos2 = chrpos2.match(/(chr\w+):(\w+)/)[2]
    if(!VV[chr1]){
            VV[chr1] = await connectVCF(chr1)
        }
        if(!VV[chr2]){ // in the rare instance chrSelect2 is different from 1
            VV[chr2] = await connectVCF(chr2)
        }
        let V1 = VV[chr1]
        let V2 = VV[chr2]
        // Extract chr values
        chr1 = chr1.match(/[0-9,X,Y,M,T]+/)[0]
        chr2 = chr2.match(/[0-9,X,Y,M,T]+/)[0]
        let q1 = await V1.query(`${chr1}:${pos1}`)
        q1.q=`${chr1}:${pos1}`
        let q2 = await V2.query(`${chr2}:${pos2}`)
        q2.q=`${chr2}:${pos2}`
        async function qq(q,V){ // making sure long rows are caught
            document.V0=V
            if(q.hit.length>0){
                console.log(`${q.q} ${q.hit.length} hit`)
                // check hit is complete
                if(V.cols.length>q.hit[0].length){ // stitch missing text
                    //let  fg = await V.fetchGz([q.ii.slice(-1)[0]-V.keyGap,q.ii.slice(-1)[0]+4*V.keyGap])
                    let txt = (await V.fetchGz(q.ii[0])).txt+(await V.fetchGz(q.ii[1])).txt
                    q.hit = txt.split('\n')
                        .filter(r=>r.match(q.q.replace(':','\t')))
                        .map(r=>r.split('\t'))
                }
            }else{
                if(!q.range){
                    q.range=q.range.dt
                }
                if(!Array.isArray(q.range)){
                    q.range=q.range.dt
                }
                if(!q.range){
                    console.log('something wrong with range for this position')
                    q.range=[]
                    //debugger
                }
                //debugger
                console.log(`${q.q} no hit`)
            }
            return q
        }
        q1 = await qq(q1,V1)
        q2 = await qq(q2,V2)
        //debugger
        //prepare counts

        let data={
            cols1:V1.cols.slice(9),
            cols2:V2.cols.slice(9),
            q1:q1,
            q2:q2,
            chrpos1:q1.q,
            chrpos2:q2.q
        }
    return data
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
        Chromossome 1: <select id="chrSelect1"></select> Position 1 <input type="number" id="pos1" value=16876455>
        <br>Chromossome 2<sup>*</sup>: <select id="chrSelect2"></select>  Position 2 <input type="number" id="pos2" value=16863727>
        <br><sub>*</sub><span style="font-size:x-small">) by default the same chromossome, different chromossomes are not linked.</span>
        <p><button id="retrievePosButton">Retrieve positions</button> directly from <a href="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" target="_blank">1000genomes</a>:</p>
        <div id="resultsDiv">...</div>`
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
        div.querySelector('#resultsDiv').innerHTML=`<span style="color:maroon">Processing ... it shouldn't take more than a minute. <br>Check console if it doesn't and report error as <a href="https://github.com/episphere/ld1000/issues" target="_blank">an issue</a>. <br>(${Date()})</span>`
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

        let chrpos1 = `chr${chr1}:${pos1}`
        let chrpos2 = `chr${chr2}:${pos2}`
        let data = await LDextraction(chrpos1,chrpos2)

        // chrpos1='chr7:16876630',chrpos2='chr7:16863828'
        
        /*
        let q1 = await V1.query(`${chr1}:${pos1}`)
        q1.q=`${chr1}:${pos1}`
        let q2 = await V2.query(`${chr2}:${pos2}`)
        q2.q=`${chr2}:${pos2}`
        async function qq(q,V){ // making sure long rows are caught
            document.V0=V
            if(q.hit.length>0){
                console.log(`${q.q} ${q.hit.length} hit`)
                // check hit is complete
                if(V.cols.length>q.hit[0].length){ // stitch missing text
                    //let  fg = await V.fetchGz([q.ii.slice(-1)[0]-V.keyGap,q.ii.slice(-1)[0]+4*V.keyGap])
                    let txt = (await V.fetchGz(q.ii[0])).txt+(await V.fetchGz(q.ii[1])).txt
                    q.hit = txt.split('\n')
                        .filter(r=>r.match(q.q.replace(':','\t')))
                        .map(r=>r.split('\t'))
                }
            }else{
                console.log(`${q.q} no hit`)
            }
            return q
        }
        q1 = await qq(q1,V1)
        q2 = await qq(q2,V2)
        //debugger
        //prepare counts

        let data={
            cols1:V1.cols.slice(9),
            cols2:V2.cols.slice(9),
            q1:q1,
            q2:q2,
            chrpos1:q1.q,
            chrpos2:q2.q
        }
        */
        
        let h2=`<hr><table><tr style="vertical-align:top"><td>`
        if(data.q1.hit.length>0){
            h2+=`1) <b style="color:maroon;font-size:large"> ${data.chrpos1}</b><span style="font-size:x-small;color:darkgreen">
            <br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp [ID]REF>ALT(QUAL-FILTER)
            <br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp [${data.q1.hit[0][2]}]${data.q1.hit[0][3]}>${data.q1.hit[0][4]}(${data.q1.hit[0][5]}-${data.q1.hit[0][6]})</span>
            <br>Close Neighbours:
            <br>${data.q1.range.map(r=>`&nbsp&nbsp&nbsp&nbsp<span style="color:blue;cursor:pointer" class="setPos1">${r[0]}:${r[1]}</span>`).join('<br>')}`
        }else{
            //debugger
            h2+=`1) <b style="color:maroon;font-size:large"> ${data.chrpos1} --> No hit !</b>
            <span style="font-size:small"><br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Maybe check kneighbouring positions.</span>
            <br>Close Neighbours:
            <br>${data.q1.range.map(r=>`&nbsp&nbsp&nbsp&nbsp<span style="color:blue;cursor:pointer" class="setPos1">${r[0]}:${r[1]}</span>`).join('<br>')}`
            //debugger
        }
        h2+=`</td><td>`
        if(data.q2.hit.length>0){
            h2+=`2) <b style="color:maroon;font-size:large"> ${data.chrpos2}</b><span style="font-size:x-small;color:darkgreen">
            <br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp [ID]REF>ALT(QUAL-FILTER)
            <br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp [${data.q2.hit[0][2]}]${data.q2.hit[0][3]}>${data.q2.hit[0][4]}(${data.q2.hit[0][5]}-${data.q2.hit[0][6]})</span>
            <br>Close Neighbours:
            <br>${data.q2.range.map(r=>`&nbsp&nbsp&nbsp&nbsp<span style="color:blue;cursor:pointer" class="setPos2">${r[0]}:${r[1]}</span>`).join('<br>')}`
        
        }else{
            h2+=`2) <b style="color:maroon;font-size:large"> ${data.chrpos2} --> No hit !</b>
            <span style="font-size:small"><br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Maybe check kneighbouring positions.</span>
            <br>Close Neighbours:
            <br>${data.q2.range.map(r=>`&nbsp&nbsp&nbsp&nbsp<span style="color:blue;cursor:pointer" class="setPos2">${r[0]}:${r[1]}</span>`).join('<br>')}`
        }
        h2+=`</td></tr><table><tr><td>`

        h2+='Counts'

        h2+=`<table id="countCombinations" style="color:maroon">`
        h2+=`<tr><td>(0|0)(0|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|0)(0|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|1)(0|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|1)(0|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|0)(1|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|0)(1|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|1)(1|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|1)(1|0)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|0)(0|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|0)(0|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|1)(0|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|1)(0|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|0)(1|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|0)(1|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(0|1)(1|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`<tr><td>(1|1)(1|1)</td><td style="color:blue">...</td><td style="font-size:x-small">...</td></tr>`
        h2+=`</table>`
        h2+=`<sub style="font-size:medium">*</sub>) <span id="OO" style="font-size:x-small">...</span>`
        
        div.querySelector('#resultsDiv')
        div.querySelector('#resultsDiv').innerHTML=h2
        let retr = div.querySelector('#retrievePosButton')
        div.querySelectorAll('.setPos1').forEach(x=>{
            let pos1 = div.querySelector('#pos1')
            x.onmouseover=function(){
                x.style.backgroundColor='yellow'
                pos1.style.backgroundColor='yellow'
            }
            x.onmouseleave=function(){
                x.style.backgroundColor=''
                pos1.style.backgroundColor=''
            }
            x.onclick=function(){
                pos1.value=this.textContent.split(':')[1]
                pos1.style.backgroundColor='yellow'
                pos1.style.color='red'
                retr.style.backgroundColor='yellow'
                retr.style.color='red'
                retr.click()
                setTimeout(function(){
                    pos1.style.backgroundColor=''
                    pos1.style.color='black'
                    retr.style.color=''
                    retr.style.backgroundColor=''
                },1000)
            }
        })
        div.querySelectorAll('.setPos2').forEach(x=>{
            let pos2 = div.querySelector('#pos2')
            x.onmouseover=function(){
                x.style.backgroundColor='yellow'
                pos2.style.backgroundColor='yellow'
            }
            x.onmouseleave=function(){
                x.style.backgroundColor=''
                pos2.style.backgroundColor=''
            }
            x.onclick=function(){
                pos2.value=this.textContent.split(':')[1]
                pos2.style.backgroundColor='yellow'
                pos2.style.color='red'
                retr.style.backgroundColor='yellow'
                retr.style.color='red'
                retr.click()
                setTimeout(function(){
                    pos2.style.backgroundColor=''
                    pos2.style.color='black'
                    retr.style.color=''
                    retr.style.backgroundColor=''
                },1000)
            }
        })

        // counting combinations if both hit

            

        if((data.q1.hit.length>0)&(data.q2.hit.length>0)){ // if both hits
            let qqPat1=data.q1.hit[0].slice(9)
            let qqPat2=data.q2.hit[0].slice(9)
            let qqPat=qqPat1.map((x,i)=>`(${x})(${qqPat2[i]})`) // combined pattern
            let tbCount = div.querySelector('#countCombinations')
            for(let i=0;i<tbCount.children[0].children.length;i++){
                let tr=tbCount.children[0].children[i]
                let pat=tr.children[0].textContent // combination pattern
                tr.children[1].textContent=qqPat.filter(x=>x==pat).length
                let pts = [] // participants
                qqPat.forEach((x,i)=>{
                    if(qqPat[i]==pat){
                        pts.push(data.cols1[i])
                    }
                })
                //if(pts.length>1000){
                //    tr.children[2].textContent='*'
                //    div.querySelector('#OO').textContent=pts.join(', ')
                //}else{
                tr.children[2].textContent=pts.slice(0,3).join(', ')
                if(pts.length>3){
                    tr.children[2].textContent+='...'
                }
                tr.children[2].onclick=function(x){
                    div.querySelector('#OO').textContent=pts.join(', ')
                }
                tr.children[2].onmouseover=function(x){
                    if(pts.length>0){
                        tr.children[2].style.backgroundColor='yellow'
                        tr.children[2].style.cursor='pointer'
                    }
                }
                tr.children[2].onmouseleave=function(x){
                    if(pts.length>0){
                        tr.children[2].style.backgroundColor=''
                        tr.children[2].style.cursor=''
                    }
                }
                //}
                //debugger
                //4
            }
        }
            
        
        //div.querySelector('#resultsDiv').innerHTML=`<hr>${JSON.stringify(q1.hit)}<hr>${JSON.stringify(q2.hit)}<hr>`
        // try 7:16876630 vs 7:16876630
        4
        

        
    }

    div.querySelector('#retrievePosButton').onclick=retrieveChrPos

    
    //retrieveChrPos()
}

export{
    connectVCF,
    UI,
    VV,
    LDextraction
}