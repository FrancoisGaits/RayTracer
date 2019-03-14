#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <iostream>
#include <stdio.h>

#include <glm/gtc/epsilon.hpp>

/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//  or boucing (add this amount to the position before casting a new ray !

const float acne_eps = 1e-4;
//const float sect_eps = 2.5e-6;

const float reflect_decay = 2.75f; //inverted evolution of clamping for reflects
const float decay_base = 3.f;

const float clamp_low = 0.f;
const float clamp_high = 7.f;

const int nb_reb = 4;

const float alias = 2;

bool intersectCylinder(Ray *ray, Intersection *intersection, Object *obj) {
  float t0, t1;

  const float cosX = dot(obj->geom.cylinder.dir,vec3(1.f,0.f,0.f)); //0
  const float cosY = dot(obj->geom.cylinder.dir,vec3(0.f,1.f,0.f)); //1
  const float cosZ = dot(obj->geom.cylinder.dir,vec3(0.f,0.f,1.f)); //0
  
  const float sinX = sqrt(1.f-cosX*cosX); //1
  const float sinY = sqrt(1.f-cosY*cosY); //0
  const float sinZ = sqrt(1.f-cosZ*cosZ); //1

  // printf("sinX : %f | sinY : %f | sinZ : %f\n",sinX,sinY,sinZ);
  
  const float cx = obj->geom.cylinder.center.x;
  const float cy = obj->geom.cylinder.center.y;
  const float cz = obj->geom.cylinder.center.z;
  
  const float dx = ray->dir.x;// * sinX;
  const float dy = ray->dir.y;// * sinY;
  const float dz = ray->dir.z;// * sinZ;

  //offsetting the ray to account the cylinder placement
  const float ox = (ray->orig.x - cx);// * sinX;
  const float oy = (ray->orig.y - cy);// * sinY;
  const float oz = (ray->orig.z - cz);// * sinZ;
  
  const float a = dx * dx * sinX + dz * dz * sinZ + (dy * dy * sinY);
  const float b = 2.f * ( dx * ox * sinX + dz  * oz * sinZ + (dy * oy * sinY));
  const float c = ox * ox * sinX + oz * oz * sinZ + (oy * oy * sinY) - obj->geom.cylinder.radius;

  const float delta = b*b - 4.f*a*c;

  if(delta < 0.f) {
    return false;
  }

  float sqrDelta = sqrt(delta);

  t0 = (-b + sqrDelta) / (2.f * a);
  t1 = (-b - sqrDelta) / (2.f * a);

  if(t0>t1) {
    intersection->inside = true;
    std::swap(t0,t1);
  }
  float y0 =  ray->orig.y * cosY + t0 * ray->dir.y * cosY + ray->orig.x*cosX + t0 * ray->dir.x*cosX + ray->orig.z*cosZ + t0 * ray->dir.z*cosZ ;
  float y1 =  ray->orig.y * cosY + t1 * ray->dir.y * cosY + ray->orig.x*cosX + t1 * ray->dir.x*cosX + ray->orig.z*cosZ + t1 * ray->dir.z*cosZ ;

  float lowerBound = obj->geom.cylinder.center.y;
  float upperBound = obj->geom.cylinder.length + lowerBound;

  if(y0<lowerBound) {
    if (y1<lowerBound) {
      return false;
    } else {
      float tPlan = t0 + (t1-t0) *(y0+upperBound) / (y0-y1);
      if (tPlan < ray->tmin || tPlan > ray->tmax) {
        return false;
      }
      ray->tmax = tPlan;
      intersection->position = rayAt(*ray,tPlan);
      intersection->mat = &obj->mat;
      intersection->normal = -obj->geom.cylinder.dir;

      return true;
    }
  } else if (y0>=lowerBound && y0<=upperBound) {
    if (t0 < ray->tmin || t0 > ray->tmax) {
      return false;
    }
   
    ray->tmax = t0;
    intersection->position = rayAt(*ray,t0);
    intersection->mat = &obj->mat;
    vec3 c = intersection->position-obj->geom.cylinder.center;
    intersection->normal = normalize(intersection->position -
				     ( (glm::dot(obj->geom.cylinder.dir,c)*obj->geom.cylinder.length)/glm::length(c) )*obj->geom.cylinder.dir);

    return true;
  } else if (y0>upperBound) {
    if(y1>upperBound) {
      return false;
    }
    else {
      float tPlan = t0 + (t1-t0) *(y0-upperBound) / (y0-y1);
      if (tPlan < ray->tmin || tPlan > ray->tmax) {
        return false;
      }
      ray->tmax = tPlan;
      intersection->position = rayAt(*ray,tPlan);
      intersection->mat = &obj->mat;
      intersection->normal = obj->geom.cylinder.dir;

      return true;
    }
  }
  return false;
}
bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
  if(glm::dot(ray->dir,obj->geom.plane.normal) == 0.f) return false;

  float t = - ( (glm::dot(ray->orig,obj->geom.plane.normal) + obj->geom.plane.dist) / (glm::dot(ray->dir,obj->geom.plane.normal)));


  if(t<ray->tmin || t>ray->tmax) return false;


  ray->tmax = t;
  intersection->position = rayAt(*ray,t);
  intersection->mat = &obj->mat;
  intersection->normal = obj->geom.plane.normal;

  return true;

}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
  vec3 dist = (ray->orig - obj->geom.sphere.center);
  /* float a = glm::dot(ray->dir,ray->dir); //si pas normalise */
  float a = 1.f;
  float b = glm::dot(ray->dir,dist) * 2.f;
  float c = glm::dot(dist,dist) - (obj->geom.sphere.radius *  obj->geom.sphere.radius) ;

  float delta = b*b - (4*a*c);
  
  if (delta < 0) return false;

  float t;

  //if(delta>-sect_eps && delta<sect_eps) {
    //return false;
  if (delta == 0) {
    t = -(0.5f*b)/a;
  } else {
    t = (-b-sqrt(delta)) / (2.f*a);
    if(t<0) {
      t = (-b+sqrt(delta)) / (2.f*a);
      intersection->inside = true;
      if(t<0) {
	return false;
      }
    }
  }
  if(t<ray->tmin || t>ray->tmax) return false;
    
  ray->tmax = t;
  intersection->position = rayAt(*ray,t);
  intersection->mat = &obj->mat;
  intersection->normal = glm::normalize(intersection->position - obj->geom.sphere.center);
 
  return true;
}

// from scratchapixel
bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj) {
  vec3 a = obj->geom.triangle.a;
  vec3 b = obj->geom.triangle.b;
  vec3 c = obj->geom.triangle.c;

  vec3 ab = b-a;
  vec3 ac = c-a;

  vec3 pvec = glm::cross(ray->dir,ac);
  float det = glm::dot(pvec,ab);

  if(det == 0.f) return false;

  float invdet = 1.f/det;
  
  vec3 tvec = ray->orig - a;

  float u = glm::dot(pvec,tvec)*invdet;
  if (u < 0 || u > 1) return false;

  vec3 qvec = glm::cross(tvec,ab);

  float v = glm::dot(ray->dir,qvec)*invdet;
  if (v < 0 || u + v > 1) return false;

  float t = glm::dot(ac,qvec)*invdet;

  if(t<ray->tmin || t>ray->tmax) return false;

  ray->tmax = t;
  intersection->position = rayAt(*ray,t);
  intersection->mat = &obj->mat;
  intersection->normal = obj->geom.triangle.normal;

  return true;
}


    
bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {

  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();
  unsigned i=0;
  Object *obj;
  while(i<objectCount) {
    obj = scene->objects[i];
    switch(obj->geom.type) {
    case SPHERE :
      hasIntersection |= intersectSphere(ray,intersection,obj);
      break;
    case PLANE :
      hasIntersection |= intersectPlane(ray,intersection,obj);
      break;
    case TRIANGLE :
      hasIntersection |= intersectTriangle(ray,intersection,obj);
      break;
    case  CYLINDER:
      hasIntersection |= intersectCylinder(ray,intersection,obj);
      break;
    default :
      hasIntersection = false;
    }

    ++i;
  }
  
  return hasIntersection;
}

/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

//limits a color component to limDown <= x <= limUp
color3 &limit(color3 &c, float limDown, float limUp) {
  c.x = c.x > limUp ? limUp : c.x;
  c.y = c.y > limUp ? limUp : c.y;
  c.z = c.z > limUp ? limUp : c.z;

  c.x = c.x < limDown ? limDown : c.x;
  c.y = c.y < limDown ? limDown : c.y;
  c.z = c.z < limDown ? limDown : c.z;

  return c;
}

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF
float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {
  float cos2 = NdotH * NdotH;
  float alpha2 = alpha * alpha;
  float tan2 = (1.f-cos2)/cos2;
  
  return RDM_chiplus(NdotH) * (std::exp((-tan2)/(alpha2))/(glm::pi<float>()*alpha2*cos2*cos2));
}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {
  float sin2tt = (pow((extIOR/intIOR),2.f) * (1.f-LdotH*LdotH));
  if (sin2tt > 1.f)
    return 1.f;

  float costt = sqrt(1.f - sin2tt);
  float Rssqrt =  ( (extIOR*LdotH) - (intIOR*costt) ) / ( (extIOR*LdotH) + (intIOR*costt) );
  float Rpsqrt =  ( (extIOR*costt) - (intIOR*LdotH) ) / ( (extIOR*costt) + (intIOR*LdotH) );
  
  return 0.5f * (Rssqrt*Rssqrt + Rpsqrt*Rpsqrt);
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {
  float tan = sqrt(1.f-DdotN*DdotN)/DdotN;
  float b = 1.f/(alpha*tan);
  float k = DdotH/DdotN;
  
  if(b>1.6f) {
    return RDM_chiplus(k);
  }

  
  return RDM_chiplus(k) * ( (3.535f*b + 2.181f*(b*b)) / (1.f + 2.276f*b + 2.577f*(b*b)) );
}


// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha) {
  
  return RDM_G1(VdotH,VdotN,alpha)*RDM_G1(LdotH,LdotN,alpha);

}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m) {
  float D = RDM_Beckmann(NdotH,m->roughness);
  float F = RDM_Fresnel(LdotH,1.f,m->IOR);
  float G = RDM_Smith(LdotH,LdotN,VdotH,VdotN,m->roughness);
  
  return m->specularColor * ( (D*F*G)/(4.f*LdotN*VdotN));

}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {
  return m->diffuseColor/glm::pi<float>();
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

  return RDM_bsdf_d(m) + RDM_bsdf_s(LdotH,NdotH,VdotH,LdotN,VdotN,m);
}



color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat) {
  color3 ret = color3{0.f};

  vec3 h = glm::normalize(v+l);
  
  float LdotN = glm::dot(l,n);
  // float LdotH = glm::dot(l,h);
  // float NdotH = glm::dot(n,h);
  // float VdotH = glm::dot(v,h);
  // float VdotN = glm::dot(v,n);
  ret = lc * RDM_bsdf(glm::dot(l,h),glm::dot(n,h),glm::dot(v,h),LdotN,glm::dot(v,n),mat) * LdotN ; 
  //ret = lc * RDM_bsdf(LdotH,NdotH,VdotH,LdotN,VdotN,mat) * LdotN; 

  return limit(ret,clamp_low,clamp_high);
}

  
//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scenem

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree) {

  color3 ret = color3{0.f};
  Intersection intersection;
  
  if(ray->depth >= nb_reb) return ret;
  
  if(intersectScene(scene,ray,&intersection)) {
    
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Lights ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(const auto &light : scene->lights) {

      vec3 l = light->position - intersection.position;
      float dist = glm::length(l);
      l /= dist;
      
      Ray shaderay;
      rayInit(&shaderay,intersection.position + acne_eps*l,l,0,dist);
      Intersection dummyinter;
      if(!intersectScene(scene,&shaderay,&dummyinter)) {
      	ret += shade(intersection.normal,-ray->dir,l,light->color,intersection.mat); //clamped by shade
      } else {
        ret += scene->ambiantLight*(intersection.mat->diffuseColor/glm::pi<float>());
      }
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Reflect ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ray reflect;
    const vec3 dir = glm::normalize(glm::reflect(ray->dir,intersection.normal));
    
    rayInit(&reflect,intersection.position + acne_eps*dir*1.15f,dir); //works with 1.15f + ??
    reflect.depth = ray->depth + 1;
    reflect.refcont = ray->refcont*reflect_decay;
    
    const vec3 h = glm::normalize(ray->dir - dir);

    float fresnel = RDM_Fresnel(glm::dot(ray->dir,h),ray->currentIOR,intersection.mat->IOR);
    
    color3 ref =  fresnel
                  * trace_ray(scene,&reflect,tree)
                  * intersection.mat->specularColor;
    if(!intersection.mat->transp)
      ret += limit(ref,clamp_low,clamp_high/(reflect.refcont+decay_base));
    // ret += RDM_Fresnel(1.f,1.f,intersection.mat->IOR) * trace_ray(scene,&reflect,tree);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Refraction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(intersection.mat->transp && fresnel < 1.f) {
      float tref = 1.f - fresnel;
      float eta;
      vec3 dirRefract;
      if(intersection.inside) {
	eta = intersection.mat->IOR; //etaMat/etaAir = etaMat
	ray->currentIOR = 1.f;
       
	dirRefract = glm::refract(ray->dir,-intersection.normal,eta);
      } else {
	eta = 1.f/intersection.mat->IOR;
	ray->currentIOR = intersection.mat->IOR;

	dirRefract = glm::refract(ray->dir,intersection.normal,eta);
      }
      Ray refractRay;
      rayInit(&refractRay,intersection.position + acne_eps*dirRefract,dirRefract,0.1);
      refractRay.depth = ray->depth + 1;

      color3 refract = tref * trace_ray(scene,&refractRay,tree);
      
      ret += limit(refract,clamp_low,clamp_high);
      
    }
    
  } else {
    ret = scene->skyColor;
  }
  
  return limit(ret,clamp_low,clamp_high);
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  float aspect = 1.f / scene->cam.aspect;

  KdTree *tree = NULL;


//! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                     aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x =
      (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;

  for (unsigned j = 0; j < img->height; ++j) {

    //~~~~~~~~~~~~~~~~~~~Affichage~~~~~~~~~~~~~~~~~~~
    if (j != 0) {
      printf("\033[A\r");
    }
    
    float progress = (float)j / img->height * 100.f;
    printf("progress\t[");
    int cpt = 0;

    for (cpt = 0; cpt < progress; cpt += 5) {
      printf(".");
    }
    
    for (; cpt < 100; cpt += 5) {
      printf(" ");
    }
    
    printf("]\n");
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
#pragma omp parallel for
    for (unsigned i = 0; i < img->width; ++i) {
      color3 *ptr = getPixelPtr(img, i, j);
      color3 pix = color3{0.f};
      
      float sub_div = alias*alias;
      
      for(size_t ali = 0; ali < alias; ++ali) {
	for(size_t alj = 0; alj < alias; ++alj) {
	
	  vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
	    float(i+(ali/alias)) * dx + float(j+(alj/alias)) * dy;

	  Ray rx;
	  rayInit(&rx, scene->cam.position, normalize(ray_dir));
	  pix += trace_ray(scene, &rx, tree);
	  
	}
      }

      *ptr = (pix/sub_div);

    }
  }
}
