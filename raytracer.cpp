
#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
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

const int nb_reb = 10;

const float alias = 3;

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
  float b = glm::dot(ray->dir,dist) * 2;
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

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {

  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();
  unsigned i=0;
  Object *obj;
  while(i<objectCount) {
    obj = scene->objects[i];
    switch(scene->objects[i]->geom.type) {
    case SPHERE :
      hasIntersection = hasIntersection | intersectSphere(ray,intersection,obj);
      break;
    case PLANE :
      hasIntersection = hasIntersection | intersectPlane(ray,intersection,obj);
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
  float LdotH = glm::dot(l,h);
  float NdotH = glm::dot(n,h);
  float VdotH = glm::dot(v,h);
  float VdotN = glm::dot(v,n);
  //  ret = lc * RDM_bsdf(glm::dot(l,h),glm::dot(n,h),glm::dot(v,h),LdotN,glm::dot(v,n),mat) * LdotN ; 
  ret = lc * RDM_bsdf(LdotH,NdotH,VdotH,LdotN,VdotN,mat) * LdotN; 
//! \todo compute bsdf, return the shaded color taking into account the
//! lightcolor


  return limit(ret,clamp_low,clamp_high);
}

  
//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scenem

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree) {

  color3 ret = color3{0.f};
  Intersection intersection;
  vec3 l;

  if(ray->depth >= nb_reb) return ret;

  
  if(intersectScene(scene,ray,&intersection)) {
    
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Lights ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(const auto &light : scene->lights) {

      l = light->position - intersection.position;
      float dist = glm::length(l);
      l /= dist;
      
      Ray shaderay;
      rayInit(&shaderay,intersection.position + acne_eps*l,l,0,dist);

      Intersection dummyinter;
      if(!intersectScene(scene,&shaderay,&dummyinter)) {
      	ret += shade(intersection.normal,-ray->dir,l,light->color,intersection.mat); //clamped by shade
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

    color3 ref =  RDM_Fresnel(glm::dot(ray->dir,h),1.f,intersection.mat->IOR)
                  * trace_ray(scene,&reflect,tree)
                  * intersection.mat->specularColor;
    
    ret += limit(ref,clamp_low,clamp_high/(reflect.refcont+decay_base));
    // ret += RDM_Fresnel(1.f,1.f,intersection.mat->IOR) * trace_ray(scene,&reflect,tree);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

  for (size_t j = 0; j < img->height; ++j) {

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
    for (size_t i = 0; i < img->width; ++i) {
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
