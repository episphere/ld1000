initStudyCase = () => {
    var chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
    var htmls=''
    chroms.forEach( e => { htmls+=`<option value="${e}" >${e}</option>` })
    ld_chrom.innerHTML=htmls
}

initStudyCase()

execute_ld = () => {
    var chrom = ld_chrom.value
    var start = ld_start.value
    var end = ld_end.value
    
    var ld_result = ld.perform_ld(chrom, start, end)
    
}
