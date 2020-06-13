afunction <- function(stuff){
  print(stuff)
}

bfunction <- function(stuff, morestuff){
  avalue <- stuff + morestuff
  return(avalue)
}

mainfunction <- function( one, two, three ){
  afunction(one)
  added <- bfunction(two,three)
  return(added)
}