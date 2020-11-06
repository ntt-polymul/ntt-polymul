print_vector(v,l) = {
  for(i=0,ceil(l/8)-1,
    printf(" ");
    for(j=1,8,printf("%7d,",v[8*i+j]));
    printf("\n");
  );
};

precomp() = {
  q = 7681;
  z = Mod(13,q);
  \\q = 10753;
  \\z = Mod(26,q);
  qinv = Mod(q^-1,2^16);
  mont = Mod(2^16,q);

  e = [384, 192, 576, 96, 480, 288, 672, 48];
  zetas = vector(16,k,centerlift(z^e[(k+1)\2]*mont));
  zetas_qinv = centerlift(qinv*zetas);
  printf("\n#define _ZETAS_PINV 128\n");
  print_vector(zetas_qinv,length(zetas_qinv));
  printf("\n#define _ZETAS 144\n");
  print_vector(zetas,length(zetas));

  \\twist level 3
  twist = vector(192,k,[]);
  f = [0, 768, e[1], e[1]+768, e[2], e[2]+768, e[3], e[3]+768,
       e[4], e[4]+768, e[5], e[5]+768, e[6], e[6]+768, e[7], e[7]+768] / 96;
  f = centerlift(Mod(f,1536/96));
  f[1] = -f[2];
  \\f = [-8, 8, 4, -4, 2, -6, 6, -2, 1, -7, 5, -3, 3, -5, 7, -1];
  for(i=0,15,
    for(j=0,5,
      l = 12*i+2*j;
      twist[l+1] = vector(16,k,centerlift(z^(f[i+1]*(16*j+k-1))*mont));
      twist[l+2] = twist[l+1];
      twist[l+1] = centerlift(twist[l+1]*qinv);
    );
  );
  twist = concat(twist);
  printf("\n#define _TWIST96 160\n");
  print_vector(twist,length(twist));

  \\twist level 6
  twist = vector(24,i,[]);
  f = [0, e[1], 768, e[1]+768, e[2], e[3], e[2]+768, e[3]+768] / 12;
  f = centerlift(Mod(f,1536/12));
  f[1] = -f[3];
  \\f = [-64, 32, 64, -32, 16, 48, -48, -16]; 
  for(i=0,3,
    for(j=0,2,
      l = 6*i+2*j;
      twist[l+1] = vector(8,k,centerlift(z^(f[2*i+(k+3)\4]*(4*j+(k+3)%4))*mont));
      twist[l+2] = twist[l+1];
      twist[l+1] = centerlift(twist[l+1]*qinv);
    );
  );
  twist = concat(twist);
  printf("\n#define _TWIST12 3232\n");
  print_vector(twist,length(twist));

  return(zetas);
}

precomp();
