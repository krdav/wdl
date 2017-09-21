import "sub.wdl" as sub
workflow test {
  String g = 'g'
  call sub.sub_wf as sc { input: thing = g }

  output {
    sc.ggg
  }
}
