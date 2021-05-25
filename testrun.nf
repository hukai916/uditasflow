process TESTRUN {
    echo true

    input:
      // val y
      tuple val(x), val(y), val(z)

    script:
      """
      echo $x, $z
      """
  }
