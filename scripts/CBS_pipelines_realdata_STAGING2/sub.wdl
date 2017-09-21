task task1 {
  String stuff

  command {
    echo 'dfgb'
  }
  output {
    String myo = 'lll'
  }
}
workflow sub_wf {
  String thing
  call task1{input: stuff=thing}

  output {
    String ggg = task1.myo
  }
}
