transcoord_func = function (si, transfunc = "gaussian", l = 0.5, c = 0, c2 = 0)
{
  si <- scale(si)
  if (transfunc == "gaussian") {
    out <- exp(-(si + c)^2/(2*l)) + c2
  }
  if (transfunc == "cosine") {
    out <- cos(pi*si/l + pi*c/2) + c2
  }
  return(out)
}
