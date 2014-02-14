function zofxy = zfn(x,y)
  global znoise;
  xvals = -1:.1:1;
  yvals = xvals;
  zxy   = sin(2.*pi.*xvals)' * sin(2.*pi.*yvals) + 0.2 .* znoise;
  zofxy = interp2(xvals, yvals, zxy, x, y);
