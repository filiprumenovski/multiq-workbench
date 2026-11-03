
/**
 *********************************************************************************************************************
 *
 *                              RE-CALCULATES THE VELOCITY VECTOR SFTRER EACH COLLISION
 *
 ************************************************************************************************************************
 */
Ion collision(Ion ion) {
  
  double vnet, lambda1, lambda2, lambda3;
  double moo1, moo2, moo3, norm1, norm2;
  double eps, theta, M, v_parallel, v_perp;
  double beta, rr1, rr2, vgas;
  double gnu1, gnu2, gnu3, norm3, vxprime, vyprime, vzprime;
  
  
  vgas = ion.minv; //terminal velocity
  vnet = Math.sqrt(ion.vx*ion.vx+ion.vy*ion.vy+ion.vz*ion.vz);
 
  if (vgas<vnet) { // re-calculate velocity in collision
    //recompute velocity vector after the collision
    // see CRC p 211, ed 28
 
    lambda1 = -ion.vy;
    moo1    =  ion.vx;
    gnu1    =  0;
 
 
    lambda2 = -ion.vx*ion.vz;
    moo2    = -ion.vy*ion.vz;
    gnu2    =  ion.vx*ion.vx+ion.vy*ion.vy;
 
    lambda3 = ion.vx;
    moo3    = ion.vy;
    gnu3    = ion.vz;

    norm1 = Math.sqrt(lambda1*lambda1 + moo1*moo1 + gnu1*gnu1);
    norm2 = Math.sqrt(lambda2*lambda2 + moo2*moo2 + gnu2*gnu2);
    norm3 = Math.sqrt(lambda3*lambda3 + moo3*moo3 + gnu3*gnu3);

    lambda1 = lambda1/norm1;
    moo1    = moo1/norm1; 
    gnu1    = gnu1/norm1;
 
    lambda2 = lambda2/norm2;
    moo2    = moo2/norm2; 
    gnu2    = gnu2/norm2;
 
    lambda3 = lambda3/norm3;
    moo3    = moo3/norm3; 
    gnu3    = gnu3/norm3;

    // compute scatering angle theta and vector projections into parallel and perp directions
    eps = random(0, 1);
    theta = Math.acos(2*eps-1);
    M = ion.m + ion.mbg;

    v_parallel = ion.mbg*Math.cos(theta)/M +ion.m/M;
    v_perp     = ion.mbg*Math.sin(theta)/M;
   
    if(vnet<vgas) { // what happens when the equilibrium reached is a complete guess work :/
        // this is an APPROXIMATION of diffusion process
       float thermoAngle = random(-2*PI,2*PI);
       float avgV = sqrt(8*1.380649*10000*273/PI*ion.MW*1.66); 
       float diffusionV = avgV+ random(0, sqrt(avgV) ); 
       v_parallel = diffusionV*Math.cos(thermoAngle);
       v_perp     = diffusionV*sin(thermoAngle);     
       // ANOTHER APPROXIMATION 
       //float thermoAngle = random(-2*pi,2*pi);
       //float diffusion = random((float)(-vgas*(1-ion.m/M)), (float)(vgas*(1-ion.m/M))); 
       //v_parallel = vgas*Math.cos(thermoAngle)+diffusion;
       //v_perp     = vgas*sin(thermoAngle)+ diffusion;
    }

    // generate second scatering angle
    beta=2*PI*random(0.0, 1.0);
    rr1=theta*180/PI;
    rr2=beta*180/PI;
    vzprime = v_parallel*vnet;
    vyprime = v_perp*vnet*Math.cos(beta);
    vxprime = v_perp*vnet*Math.sin(beta);

    // re-transform back to original coordinate system
    ion.vx=lambda1*vxprime+lambda2*vyprime+lambda3*vzprime;
    ion.vy=moo1*vxprime+moo2*vyprime+moo3*vzprime;
    ion.vz=gnu1*vxprime+gnu2*vyprime+gnu3*vzprime;
 
  } // close in vnet>vgas
  
  return ion;
 
}
