make_test_data <- function() {
  # 3 treat, 2 control
  treat1 <- tibble::tibble(
    chr = "chr1", pos = 1:6, ref = c("A","A","T","C","A","G"),
    A = c(100, 200, 0, 0, 100, 0),
    T = c(  0,   0, 80, 0,   0, 0),
    C = c(  0,   0, 20, 50,  0, 0),
    G = c( 15,   0,  0, 0,   3, 60),
    site_id = paste0("chr1_", 1:6)
  )
  treat2 <- tibble::tibble(
    chr = "chr1", pos = 1:6, ref = c("A","A","T","C","A","G"),
    A = c( 90, 180, 0, 0, 110, 0),
    T = c(  0,   0, 75, 0,   0, 0),
    C = c(  0,   0, 25, 55,  0, 0),
    G = c( 12,   0,  0, 0,   0, 58),
    site_id = paste0("chr1_", 1:6)
  )
  treat3 <- tibble::tibble(
    chr = "chr1", pos = 1:6, ref = c("A","A","T","C","A","G"),
    A = c( 95, 190, 0, 0, 105, 0),
    T = c(  0,   0, 70, 0,   0, 0),
    C = c(  0,   0, 30, 48,  0, 0),
    G = c( 10,   0,  0, 0,   0, 62),
    site_id = paste0("chr1_", 1:6)
  )
  cont1 <- tibble::tibble(
    chr = "chr1", pos = 1:6, ref = c("A","A","T","C","A","G"),
    A = c(120, 210, 0, 0, 115, 0),
    T = c(  0,   0, 95, 0,   0, 0),
    C = c(  0,   0,  5, 60,  0, 0),
    G = c(  2,   0,  0, 0,   0, 55),
    site_id = paste0("chr1_", 1:6)
  )
  cont2 <- tibble::tibble(
    chr = "chr1", pos = 1:6, ref = c("A","A","T","C","A","G"),
    A = c(130, 205, 0, 0, 120, 0),
    T = c(  0,   0, 90, 0,   0, 0),
    C = c(  0,   0,  8, 52,  0, 0),
    G = c(  1,   0,  0, 0,   0, 50),
    site_id = paste0("chr1_", 1:6)
  )

  list(T1 = treat1, T2 = treat2, T3 = treat3, C1 = cont1, C2 = cont2)
}


#
#
# | Site | Ref | Edit | Treat G/C total | Reps>0 | ref!=target | Keep? |
# |------|-----|------|-----------------|--------|------------|-------|
# | 1    | A   | A>G  | 15+12+10=37     | 3      | A!=G       |   YES |
# | 2    | A   | A>G  | 0               | 0      | ✓          |   NO  |
# | 3    | T   | T>C  | 20+25+30=75     | 3      | T!=C       |   YES |
# | 4    | C   |  —   | not A or T ref  | —      | —          |   NO  |
# | 5    | A   | A>G  | 3+0+0=3.        | 1      | ✓          |   NO  |
# | 6    | G   | A>G  | —               | —      | G==G       |   NO  |
#
#

