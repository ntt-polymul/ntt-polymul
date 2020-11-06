brv = [128, 64, 192, 32, 160, 96, 224, 16, 144, 80, 208, 48, 176, 112, 240, 8, 136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248, 4, 132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244, 12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 188, 124, 252, 2, 130, 66, 194, 34, 162, 98, 226, 18, 146, 82, 210, 50, 178, 114, 242, 10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250, 6, 134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246, 14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94, 222, 62, 190, 126, 254, 1, 129, 65, 193, 33, 161, 97, 225, 17, 145, 81, 209, 49, 177, 113, 241, 9, 137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249, 5, 133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245, 13, 141, 77, 205, 45, 173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253, 3, 131, 67, 195, 35, 163, 99, 227, 19, 147, 83, 211, 51, 179, 115, 243, 11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251, 7, 135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247, 15, 143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255];

print_vector(v,l) = {
  for(i=0,ceil(l/8)-1,
    printf(" ");
    for(j=1,8,printf("%7d,",v[8*i+j]));
    printf("\n");
  );
};

precomp() = {
  \\q = 3329;
  \\z = Mod(17,q);
  \\q = 7681;
  \\z = Mod(62,q);
  q = 10753;
  z = Mod(10,q);
  qinv = Mod(q^-1,2^16);
  mont = Mod(2^16,q);

  zetas_exp = vector(10, i, []);

  \\ level 0
  zetas_exp[2] = vector(16, i, centerlift(z^brv[1]*mont));
  zetas_exp[1] = centerlift(qinv*zetas_exp[2]);

  \\levels 1-2
  k = 3;
  for(level = 1,2,
    rep = 16/2^(level-1);
    for(i = 0,2^(level-1)-1,
      zetas_exp[k] = concat(zetas_exp[k],vector(rep,k,centerlift(z^brv[2^level+i]*mont)));
    );
    zetas_exp[k+1] = zetas_exp[k+0];
    zetas_exp[k+0] = centerlift(zetas_exp[k+0]*qinv);
    k += 2;
  );
  for(level = 1,2,
    rep = 16/2^(level-1);
    for(i = 0,2^(level-1)-1,
      zetas_exp[k] = concat(zetas_exp[k],vector(rep,k,centerlift(z^brv[3*2^(level-1)+i]*mont)));
    );
    zetas_exp[k+1] = zetas_exp[k+0];
    zetas_exp[k+0] = centerlift(zetas_exp[k+0]*qinv);
    k += 2;
  );
  zetas_exp = concat(zetas_exp);
  printf("\n#define _ZETAS 128\n");
  print_vector(zetas_exp,length(zetas_exp));

  \\twist level 3
  twist = vector(32,i,[]);
  e = vector(8,k,centerlift(Mod(brv[8+k-1]/16,512/32)));
  \\ e = [1, -7, 5, -3, 3, -5, 7, -1];
  for(i=0,7,
    for(j=1,4,
      twist[2*i+1] = concat(twist[2*i+1],vector(4,k,centerlift(z^(e[j]*(4*i+k-1))*mont)));
    );
    twist[2*i+2] = twist[2*i+1];
    twist[2*i+1] = centerlift(twist[2*i+1]*qinv);
  );
  for(i=0,7,
    for(j=5,8,
      twist[16+2*i+1] = concat(twist[16+2*i+1], vector(4,k,centerlift(z^(e[j]*(4*i+k-1))*mont)));
    );
    twist[16+2*i+2] = twist[16+2*i+1];
    twist[16+2*i+1] = centerlift(twist[16+2*i+1]*qinv);
  );
  twist = concat(twist);
  printf("\n#define _TWIST32 288\n");
  print_vector(twist,length(twist));

  \\twist level 6
  e = [0, 256, brv[1], brv[1]+256, brv[2], brv[2]+256, brv[3], brv[3]+256] / 4;
  e = centerlift(Mod(e,512/4));
  e[1] = -e[2];
  \\ e = [-64, 64, 32, -32, 16, -48, 48, -16];
  twist = vector(16,i,[]);
  twist[ 1] = vector(4,k,centerlift(z^(e[1]*(k-1))*mont));
  twist[ 3] = vector(4,k,centerlift(z^(e[2]*(k-1))*mont));
  twist[ 5] = vector(4,k,centerlift(z^(e[3]*(k-1))*mont));
  twist[ 7] = vector(4,k,centerlift(z^(e[4]*(k-1))*mont));
  twist[ 9] = vector(4,k,centerlift(z^(e[5]*(k-1))*mont));
  twist[11] = vector(4,k,centerlift(z^(e[6]*(k-1))*mont));
  twist[13] = vector(4,k,centerlift(z^(e[7]*(k-1))*mont));
  twist[15] = vector(4,k,centerlift(z^(e[8]*(k-1))*mont));
  twist[2] = twist[1];
  twist[1] = centerlift(twist[1]*qinv);
  twist[4] = twist[3];
  twist[3] = centerlift(twist[3]*qinv);
  twist[6] = twist[5];
  twist[5] = centerlift(twist[5]*qinv);
  twist[8] = twist[7];
  twist[7] = centerlift(twist[7]*qinv);
  twist[10] = twist[9];
  twist[9] = centerlift(twist[9]*qinv);
  twist[12] = twist[11];
  twist[11] = centerlift(twist[11]*qinv);
  twist[14] = twist[13];
  twist[13] = centerlift(twist[13]*qinv);
  twist[16] = twist[15];
  twist[15] = centerlift(twist[15]*qinv);
  twist = concat(twist);
  printf("\n#define _TWIST4 800\n");
  print_vector(twist,length(twist));
}

precomp();
