print_vector(v,l) = {
  for(i=0,ceil(l/8)-1,
    printf(" ");
    for(j=1,8,printf("%7d,",v[8*i+j]));
    printf("\n");
  );
};

precomp() = {
  \\q = 3457;
  \\z = Mod(3,q);
  q = 8641;
  z = Mod(14,q);
  qinv = Mod(q^-1,2^16);
  mont = Mod(2^16,q);

  e = [576,192,384,768];
  zetas3 = vector(8,k,centerlift(z^e[(k+1)\2]*mont));
  zetas3_qinv = centerlift(qinv*zetas3);
  printf("\n#define _ZETAS3_PINV 128\n");
  print_vector(zetas3_qinv,length(zetas3_qinv));
  printf("\n#define _ZETAS3 136\n");
  print_vector(zetas3,length(zetas3));
  zetas3_inv = vector(8,k,centerlift(z^-e[(k+1)\2]*mont));
  zetas3_inv_qinv = centerlift(qinv*zetas3_inv);
  printf("\n#define _ZETAS3_INV_PINV 144\n");
  print_vector(zetas3_inv_qinv,length(zetas3_inv_qinv));
  printf("\n#define _ZETAS3_INV 152\n");
  print_vector(zetas3_inv,length(zetas3_inv));

  f = [432,216,648,108];
  zetas = vector(8,k,centerlift(z^f[(k+1)\2]*mont));
  zetas_qinv = centerlift(qinv*zetas);
  printf("\n#define _ZETAS_PINV 160\n");
  print_vector(zetas_qinv,length(zetas_qinv));
  printf("\n#define _ZETAS 168\n");
  print_vector(zetas,length(zetas));

  \\twist level 1
  twist = vector(216,i,[]);
  g = [0, e[1], 2*e[1], e[2], e[2]+e[1], e[2]+2*e[1], e[3], e[3]+e[1], e[3]+2*e[1]] / 192;
  g = centerlift(Mod(g,1728/192));
  \\g = [0, 3, -3, 1, 4, -2, 2, -4, -1];
  for(j=0,11,
    for(i=0,8,
      l = 18*j+2*i;
      twist[l+1] = vector(16,k,centerlift(z^(g[i+1]*(16*j+k-1))*mont));
      twist[l+2] = twist[l+1];
      twist[l+1] = centerlift(twist[l+1] * qinv);
    );
  );
  twist = concat(twist);
  printf("\n#define _TWIST192 176\n");
  print_vector(twist,length(twist));

  \\twist level 5
  twist = vector(24,i,[]);
  g = [0, 864, f[1], f[1]+864, f[2], f[2]+864, f[3], f[3]+864] / 24;
  g = centerlift(Mod(g,1728/24));
  twist[ 1] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 0))*mont));
  twist[ 3] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 2))*mont));
  twist[ 5] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 4))*mont));
  twist[ 7] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 6))*mont));
  twist[ 9] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 8))*mont));
  twist[11] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+10))*mont));
  twist[13] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+12))*mont));
  twist[15] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+14))*mont));
  twist[17] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+16))*mont));
  twist[19] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+18))*mont));
  twist[21] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+20))*mont));
  twist[23] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+22))*mont));
  for(i=0,11,
    twist[2*i+2] = twist[2*i+1];
    twist[2*i+1] = centerlift(twist[2*i+1]*qinv);
  );
  twist = concat(twist);
  printf("\n#define _TWIST24 3632\n");
  print_vector(twist,length(twist));

  g = -g;
  twist[ 1] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 0))*mont));
  twist[ 3] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 2))*mont));
  twist[ 5] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 4))*mont));
  twist[ 7] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 6))*mont));
  twist[ 9] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+ 8))*mont));
  twist[11] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+10))*mont));
  twist[13] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+12))*mont));
  twist[15] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+14))*mont));
  twist[17] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+16))*mont));
  twist[19] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+18))*mont));
  twist[21] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+20))*mont));
  twist[23] = vector(16,k,centerlift(z^(g[(k+1)\2]*((k+1)%2+22))*mont));
  for(i=0,11,
    twist[2*i+2] = twist[2*i+1];
    twist[2*i+1] = centerlift(twist[2*i+1]*qinv);
  );
  twist = concat(twist);
  printf("\n#define _TWIST24INV 4016\n");
  print_vector(twist,length(twist));

  return(zetas);
}

precomp();
