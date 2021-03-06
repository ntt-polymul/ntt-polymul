brv = [0, 128, 64, 192, 32, 160, 96, 224, 16, 144, 80, 208, 48, 176, 112, 240, 8, 136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248, 4, 132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244, 12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 188, 124, 252, 2, 130, 66, 194, 34, 162, 98, 226, 18, 146, 82, 210, 50, 178, 114, 242, 10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250, 6, 134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246, 14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94, 222, 62, 190, 126, 254, 1, 129, 65, 193, 33, 161, 97, 225, 17, 145, 81, 209, 49, 177, 113, 241, 9, 137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249, 5, 133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245, 13, 141, 77, 205, 45, 173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253, 3, 131, 67, 195, 35, 163, 99, 227, 19, 147, 83, 211, 51, 179, 115, 243, 11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251, 7, 135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247, 15, 143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255];

print_vector(v,l) = {
  for(i=0,ceil(l/8)-1,
    printf(" ");
    for(j=1,8,printf("%7d,",v[8*i+j]));
    printf("\n");
  );
};

precomp() = {
  \\q = 7681;
  \\z = Mod(62,q);
  q = 10753;
  z = Mod(10,q);
  qinv = Mod(q^-1,2^16);
  mont = Mod(2^16,q);

  zetas = vector(256,k,centerlift(z^brv[k]*mont));
  zetas_qinv = centerlift(qinv*zetas);
  printf("\n#define _ZETAS_PINV 128\n");
  print_vector(zetas_qinv,length(zetas_qinv));
  printf("\n#define _ZETAS 384\n");
  print_vector(zetas,length(zetas));

  zetas = vector(8,k,centerlift(z^-brv[k]*mont));
  zetas_qinv = centerlift(qinv*zetas);
  printf("\n#define _ZETAS_INV_PINV 640\n");
  print_vector(zetas_qinv,length(zetas_qinv));
  printf("\n#define _ZETAS_INV 648\n");
  print_vector(zetas,length(zetas));
}

precomp();
