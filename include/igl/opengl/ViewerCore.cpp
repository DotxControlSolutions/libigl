// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "ViewerCore.h"
#include "ViewerData.h"
#include "gl.h"
#include "../quat_to_mat.h"
#include "../snap_to_fixed_up.h"
#include "../look_at.h"
#include "../frustum.h"
#include "../ortho.h"
#include "../massmatrix.h"
#include "../barycenter.h"
#include "../PI.h"
#include <Eigen/Geometry>
#include <iostream>

IGL_INLINE void igl::opengl::ViewerCore::align_camera_center(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F)
{
  if(V.rows() == 0)
    return;

  //get_scale_and_shift_to_fit_mesh(V, F, model_zoom, model_translation);
  get_scale_and_shift_to_fit_mesh(V,F,camera_base_zoom,camera_base_translation);
  // Rather than crash on empty mesh...
  if(V.size() > 0)
  {
    object_scale = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  }
}

IGL_INLINE void igl::opengl::ViewerCore::get_scale_and_shift_to_fit_mesh(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  float& zoom,
  Eigen::Vector3f& shift)
{
  if (V.rows() == 0)
    return;

  Eigen::MatrixXd BC;
  if (F.rows() <= 1)
  {
    BC = V;
  } else
  {
    igl::barycenter(V,F,BC);
  }
  return get_scale_and_shift_to_fit_mesh(BC,zoom,shift);
}

IGL_INLINE void igl::opengl::ViewerCore::align_camera_center(
  const Eigen::MatrixXd& V)
{
  if(V.rows() == 0)
    return;

  //get_scale_and_shift_to_fit_mesh(V, model_zoom, model_translation);
  get_scale_and_shift_to_fit_mesh(V,camera_base_zoom,camera_base_translation);
  // Rather than crash on empty mesh...
  if(V.size() > 0)
  {
    object_scale = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  }
}

IGL_INLINE void igl::opengl::ViewerCore::get_scale_and_shift_to_fit_mesh(
  const Eigen::MatrixXd& V,
  float& zoom,
  Eigen::Vector3f& shift)
{
  if (V.rows() == 0)
    return;

  auto min_point = V.colwise().minCoeff();
  auto max_point = V.colwise().maxCoeff();
  auto centroid  = (0.5*(min_point + max_point)).eval();
  shift.setConstant(0);
  shift.head(centroid.size()) = -centroid.cast<float>();
  zoom = 2.0 / (max_point-min_point).array().abs().maxCoeff();
}


IGL_INLINE void igl::opengl::ViewerCore::clear_framebuffers()
{
  // The glScissor call ensures we only clear this core's buffers,
  // (in case the user wants different background colors in each viewport.)
  glScissor(viewport(0), viewport(1), viewport(2), viewport(3));
  glEnable(GL_SCISSOR_TEST);
  glClearColor(background_color[0],
               background_color[1],
               background_color[2],
               background_color[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glDisable(GL_SCISSOR_TEST);
}

IGL_INLINE void igl::opengl::ViewerCore::draw(
  ViewerData& data,
  bool update_matrices)
{
  
  using namespace std;
  using namespace Eigen;
  if (depth_test)
  {
      glEnable(GL_DEPTH_TEST);
      glDepthFunc(GL_LESS);
  }
  else
    glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  
  /* Bind and potentially refresh mesh/line/point data */
  if (data.dirty)
  {
    data.updateGL(data, data.invert_normals, data.meshgl);
    data.dirty = MeshGL::DIRTY_NONE;
  }
  
  data.meshgl.bind_mesh();

  // Initialize uniform
  glViewport(viewport(0), viewport(1), viewport(2), viewport(3));

  if(update_matrices)
  {
    //model = Eigen::Matrix4f::Identity();
    view = Eigen::Matrix4f::Identity();
    proj = Eigen::Matrix4f::Identity();
    //norm = Eigen::Matrix4f::Identity();

    //std::cout << "Viewer.core().draw(), line 153" << std::endl;

    // Set view
    look_at( camera_eye, camera_center, camera_up, view);
    view = view
      * (trackball_angle * Eigen::Scaling(camera_zoom * camera_base_zoom)
      * Eigen::Translation3f(camera_translation + camera_base_translation)).matrix();
    //cout << "camera_center: " << camera_center << endl;
    //cout << "camera_eye: " << camera_eye << endl;
    //cout << "camera_zoom: " << camera_zoom << endl;
    //cout << "camera_base_zoom: " << camera_base_zoom << endl;
    //cout << "camera translation: " << camera_translation << endl; // <- this one should move when the right mouse button is used
    //cout << "camera_base_translation: " << camera_base_translation << endl;

    float width = viewport(2);
    float height = viewport(3);
    //norm = view.inverse().transpose();

    // Set projection
    if (orthographic)
    {
      float length = (camera_eye - camera_center).norm();
      float h = tan(camera_view_angle/360.0 * igl::PI) * (length);
      ortho(-h*width/height, h*width/height, -h, h, camera_dnear, camera_dfar,proj);
    }
    else
    {
      float fH = tan(camera_view_angle / 360.0 * igl::PI) * camera_dnear;
      float fW = fH * (double)width/(double)height;
      frustum(-fW, fW, -fH, fH, camera_dnear, camera_dfar,proj);
    }
    // end projection

    //// Set model transformation
    //float mat[16];
    //igl::quat_to_mat(trackball_angle.coeffs().data(), mat);

    //for (unsigned i = 0; i < 4; ++i)
    //    for (unsigned j = 0; j < 4; ++j)
    //        model(i, j) = mat[i + 4 * j];

    //// Why not just use Eigen::Transform<double,3,Projective> for model...?
    //model.topLeftCorner(3, 3) *= camera_zoom;
    //model.topLeftCorner(3, 3) *= model_zoom;
    //model.col(3).head(3) += model.topLeftCorner(3, 3) * model_translation;
  }


  // Send transformations to the GPU
  GLint modeli = glGetUniformLocation(data.meshgl.shader_mesh,"model"); // Add the model matrix
  GLint viewi  = glGetUniformLocation(data.meshgl.shader_mesh,"view");
  GLint proji  = glGetUniformLocation(data.meshgl.shader_mesh,"proj");
 // GLint normi  = glGetUniformLocation(data.meshgl.shader_mesh,"normal_matrix");

  model = data.rotation_matrix;

  //cout << "view: " << view << endl;
  //cout << "model: " << model << endl;
  //cout << "proj: " << proj << endl;

  glUniformMatrix4fv(modeli, 1, GL_FALSE, data.rotation_matrix.data());
  glUniformMatrix4fv(viewi , 1, GL_FALSE, view.data());
  glUniformMatrix4fv(proji , 1, GL_FALSE, proj.data());
  //glUniformMatrix4fv(normi , 1, GL_FALSE, norm.data());
  
  // Light parameters
  GLint specular_exponenti    = glGetUniformLocation(data.meshgl.shader_mesh,"specular_exponent");
  GLint key_light_i = glGetUniformLocation(data.meshgl.shader_mesh, "key_light");
  GLint fill_light_i = glGetUniformLocation(data.meshgl.shader_mesh, "fill_light");
  GLint back_light_i = glGetUniformLocation(data.meshgl.shader_mesh, "back_light");
  GLint light_position_worldi = glGetUniformLocation(data.meshgl.shader_mesh, "light_position_world");
  GLint light_position_world2i = glGetUniformLocation(data.meshgl.shader_mesh, "light_position_world2");
  //GLint light_position_world3i = glGetUniformLocation(data.meshgl.shader_mesh, "light_position_world3");
  //GLint light_position_world4i = glGetUniformLocation(data.meshgl.shader_mesh, "light_position_world4");
  //GLint light_position_eyei   = glGetUniformLocation(data.meshgl.shader_mesh,"light_position_eye");
  GLint lighting_factori      = glGetUniformLocation(data.meshgl.shader_mesh,"lighting_factor");
  GLint fixed_colori          = glGetUniformLocation(data.meshgl.shader_mesh,"fixed_color");
  GLint texture_factori       = glGetUniformLocation(data.meshgl.shader_mesh,"texture_factor");
  GLint matcap_factori        = glGetUniformLocation(data.meshgl.shader_mesh,"matcap_factor");
  GLint double_sidedi         = glGetUniformLocation(data.meshgl.shader_mesh,"double_sided");

  // compute camera position in world-space
  Eigen::Matrix<float, 4, 1> camera_pos_4;
  camera_pos_4 << camera_eye, 1.0f;
  Eigen::Vector3f camera_pos_ws = (view.inverse() * camera_pos_4).cast<float>().head(3);

  float distance_cam_to_center = (camera_pos_ws - camera_center).norm();

  // compute camera reference frame in world space - up
  Eigen::Matrix<float, 4, 1> camera_up_4;
  camera_up_4 << camera_up, 0.0f;
  Eigen::Vector3f camera_up_ws = (view.inverse() * camera_up_4).cast<float>().head(3);
  camera_up_ws.normalize();

  // compute camera reference frame in world space - dir
  Eigen::Matrix<float, 4, 1> camera_dir_4;
  camera_dir_4 << camera_dir, 0.0f;
  Eigen::Vector3f camera_dir_ws = (view.inverse() * camera_dir_4).cast<float>().head(3);
  camera_dir_ws.normalize();

  // compute camera reference frame in world space - side
  Eigen::Matrix<float, 4, 1> camera_side_4;
  camera_side_4 << camera_side, 0.0f;
  Eigen::Vector3f camera_side_ws = (view.inverse() * camera_side_4).cast<float>().head(3);
  camera_side_ws.normalize();

  // direction camera to back light
  Eigen::Vector3f dir_cam_to_back_light;
  dir_cam_to_back_light = -1*camera_side_ws + 2*camera_dir_ws + 2*camera_up_ws;
  dir_cam_to_back_light.normalize();

  Eigen::Vector3f key_light = camera_pos_ws + 0.5 * distance_cam_to_center * camera_side_ws;
  Eigen::Vector3f fill_light = camera_pos_ws - 1.0 * distance_cam_to_center * camera_side_ws;
  Eigen::Vector3f back_light = camera_pos_ws + 3 * distance_cam_to_center * dir_cam_to_back_light;

  //cout << "camera_pos_ws: " << camera_pos_ws << endl;
  //cout << "distance_cam_to_center: " << distance_cam_to_center << endl;
  //cout << "key_light: " << key_light << endl;
  //cout << "fill_light: " << fill_light << endl;
  //cout << "back_light: " << back_light << endl;

  glUniform1f(specular_exponenti, data.shininess);
  /*glUniform3fv(light_position_eyei, 1, light_position.data());
  Vector3f rev_light = -1. * light_position;
  glUniform3fv(light_position_worldi, 1, rev_light.data());*/
  glUniform3fv(key_light_i, 1, key_light.data());
  glUniform3fv(fill_light_i, 1, fill_light.data());
  glUniform3fv(back_light_i, 1, back_light.data());
  glUniform3fv(light_position_worldi, 1, light_position.data());
  glUniform3fv(light_position_world2i, 1, light_position_2.data());
  //glUniform3fv(light_position_world3i, 1, light_position_3.data());
  //glUniform3fv(light_position_world4i, 1, light_position_4.data());
  glUniform1f(lighting_factori, lighting_factor); // enables lighting
  glUniform4f(fixed_colori, 0.0, 0.0, 0.0, 0.0);
  

  //std::cout << "Viewer.core().draw(), line 239" << std::endl;
  if (data.V.rows()>0)
  {
    // Render fill
    if (is_set(data.show_faces))
    {
      // Texture
      glUniform1f(texture_factori, is_set(data.show_texture) ? 1.0f : 0.0f);
      glUniform1f(matcap_factori, is_set(data.use_matcap) ? 1.0f : 0.0f);
      glUniform1f(double_sidedi, data.double_sided ? 1.0f : 0.0f);
      //std::cout << "Viewer.core().draw(), line 245" << std::endl;
      data.meshgl.draw_mesh(true);
      glUniform1f(matcap_factori, 0.0f);
      glUniform1f(texture_factori, 0.0f);
    }

    // Render wireframe
    if (is_set(data.show_lines))
    {
      glLineWidth(data.line_width);
      glUniform4f(fixed_colori,
        data.line_color[0],
        data.line_color[1],
        data.line_color[2],
        data.line_color[3]);
      //std::cout << "Viewer.core().draw(), line 260" << std::endl;
      data.meshgl.draw_mesh(false);
      glUniform4f(fixed_colori, 0.0f, 0.0f, 0.0f, 0.0f);
    }
  }

  if (is_set(data.show_overlay))
  {
    if (is_set(data.show_overlay_depth))
      glEnable(GL_DEPTH_TEST);
    else
      glDisable(GL_DEPTH_TEST);

    if (data.lines.rows() > 0)
    {
      data.meshgl.bind_overlay_lines();
      modeli = glGetUniformLocation(data.meshgl.shader_overlay_lines, "model");
      viewi  = glGetUniformLocation(data.meshgl.shader_overlay_lines,"view");
      proji  = glGetUniformLocation(data.meshgl.shader_overlay_lines,"proj");

      glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
      // This must be enabled, otherwise glLineWidth has no effect
      glEnable(GL_LINE_SMOOTH);
      glLineWidth(data.line_width);

      data.meshgl.draw_overlay_lines();
    }

    if (data.points.rows() > 0)
    {
      data.meshgl.bind_overlay_points();
      modeli = glGetUniformLocation(data.meshgl.shader_overlay_points, "model");
      viewi  = glGetUniformLocation(data.meshgl.shader_overlay_points,"view");
      proji  = glGetUniformLocation(data.meshgl.shader_overlay_points,"proj");

      glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
      glPointSize(data.point_size);
      data.meshgl.draw_overlay_points();
    }
    glEnable(GL_DEPTH_TEST);
  }

  if(is_set(data.show_vertex_labels)&&data.vertex_labels_positions.rows()>0) 
    draw_labels(data, data.meshgl.vertex_labels);
  if(is_set(data.show_face_labels)&&data.face_labels_positions.rows()>0) 
    draw_labels(data, data.meshgl.face_labels);
  if(is_set(data.show_custom_labels)&&data.labels_positions.rows()>0) 
    draw_labels(data, data.meshgl.custom_labels);
}

IGL_INLINE void igl::opengl::ViewerCore::draw_buffer(
  ViewerData& data,
  bool update_matrices,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A)
{
  assert(R.rows() == G.rows() && G.rows() == B.rows() && B.rows() == A.rows());
  assert(R.cols() == G.cols() && G.cols() == B.cols() && B.cols() == A.cols());

  unsigned width = R.rows();
  unsigned height = R.cols();
  if(width == 0 && height == 0)
  {
    width = viewport(2);
    height = viewport(3);
  }
  R.resize(width,height);
  G.resize(width,height);
  B.resize(width,height);
  A.resize(width,height);

  ////////////////////////////////////////////////////////////////////////
  // PREPARE width×height BUFFERS does *not* depend on `data`
  //   framebuffer
  //   textureColorBufferMultiSampled
  //   rbo
  //   intermediateFBO
  //   screenTexture
  //
  ////////////////////////////////////////////////////////////////////////
  // https://learnopengl.com/Advanced-OpenGL/Anti-Aliasing

  // Create an initial multisampled framebuffer
  unsigned int framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
  // create a multisampled color attachment texture
  unsigned int textureColorBufferMultiSampled;
  glGenTextures(1, &textureColorBufferMultiSampled);
  glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, textureColorBufferMultiSampled);
  glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, GL_RGBA, width, height, GL_TRUE);
  glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, textureColorBufferMultiSampled, 0);
  // create a (also multisampled) renderbuffer object for depth and stencil attachments
  unsigned int rbo;
  glGenRenderbuffers(1, &rbo);
  glBindRenderbuffer(GL_RENDERBUFFER, rbo);
  glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_DEPTH24_STENCIL8, width, height);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);
  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // configure second post-processing framebuffer
  unsigned int intermediateFBO;
  glGenFramebuffers(1, &intermediateFBO);
  glBindFramebuffer(GL_FRAMEBUFFER, intermediateFBO);
  // create a color attachment texture
  unsigned int screenTexture;
  glGenTextures(1, &screenTexture);
  glBindTexture(GL_TEXTURE_2D, screenTexture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, screenTexture, 0);	// we only need a color buffer
  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // Clear the buffer
  glClearColor(background_color(0), background_color(1), background_color(2), 0.f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // Save old viewport
  Eigen::Vector4f viewport_ori = viewport;
  viewport << 0,0,width,height;
  // Draw
  draw(data,update_matrices);
  // Restore viewport
  viewport = viewport_ori;

  glBindFramebuffer(GL_READ_FRAMEBUFFER, framebuffer);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, intermediateFBO);
  glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST);

  glBindFramebuffer(GL_FRAMEBUFFER, intermediateFBO);
  // Copy back in the given Eigen matrices
  GLubyte* pixels = (GLubyte*)calloc(width*height*4,sizeof(GLubyte));
  glReadPixels(0, 0,width, height,GL_RGBA, GL_UNSIGNED_BYTE, pixels);

  // Clean up
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
  glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glDeleteTextures(1, &screenTexture);
  glDeleteTextures(1, &textureColorBufferMultiSampled);
  glDeleteFramebuffers(1, &framebuffer);
  glDeleteFramebuffers(1, &intermediateFBO);
  glDeleteRenderbuffers(1, &rbo);

  int count = 0;
  for (unsigned j=0; j<height; ++j)
  {
    for (unsigned i=0; i<width; ++i)
    {
      R(i,j) = pixels[count*4+0];
      G(i,j) = pixels[count*4+1];
      B(i,j) = pixels[count*4+2];
      A(i,j) = pixels[count*4+3];
      ++count;
    }
  }
  // Clean up
  free(pixels);
}

// Define uniforms for text labels
IGL_INLINE void igl::opengl::ViewerCore::draw_labels(
  ViewerData& data,
  const igl::opengl::MeshGL::TextGL& labels
){
  glDisable(GL_LINE_SMOOTH); // Clear settings if overlay is activated
  data.meshgl.bind_labels(labels);
  GLint viewi = glGetUniformLocation(data.meshgl.shader_text,"view");
  GLint proji = glGetUniformLocation(data.meshgl.shader_text,"proj");
  glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
  glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
  // Parameters for mapping characters from font atlass
  float width  = viewport(2);
  float height = viewport(3);
  float text_shift_scale_factor = orthographic ? 0.01 : 0.03;
  float render_scale = (orthographic ? 0.6 : 1.7) * data.label_size;
  glUniform1f(glGetUniformLocation(data.meshgl.shader_text, "TextShiftFactor"), text_shift_scale_factor);
  glUniform3f(glGetUniformLocation(data.meshgl.shader_text, "TextColor"), data.label_color(0), data.label_color(1), data.label_color(2));
  glUniform2f(glGetUniformLocation(data.meshgl.shader_text, "CellSize"), 1.0f / 16, (300.0f / 384) / 6);
  glUniform2f(glGetUniformLocation(data.meshgl.shader_text, "CellOffset"), 0.5 / 256.0, 0.5 / 256.0);
  glUniform2f(glGetUniformLocation(data.meshgl.shader_text, "RenderSize"), 
                                    render_scale * 0.75 * 16 / (width), 
                                    render_scale * 0.75 * 33.33 / (height));
  glUniform2f(glGetUniformLocation(data.meshgl.shader_text, "RenderOrigin"), -2, 2);
  data.meshgl.draw_labels(labels);
  glEnable(GL_DEPTH_TEST);
}

IGL_INLINE void igl::opengl::ViewerCore::set_rotation_type(
  const igl::opengl::ViewerCore::RotationType & value)
{
  using namespace Eigen;
  using namespace std;
  const RotationType old_rotation_type = rotation_type;
  rotation_type = value;
  if(rotation_type == ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP &&
    old_rotation_type != ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
  {
    snap_to_fixed_up(Quaternionf(trackball_angle),trackball_angle);
  }
}

IGL_INLINE void igl::opengl::ViewerCore::set(unsigned int &property_mask, bool value) const
{
  if (!value)
    unset(property_mask);
  else
    property_mask |= id;
}

IGL_INLINE void igl::opengl::ViewerCore::unset(unsigned int &property_mask) const
{
  property_mask &= ~id;
}

IGL_INLINE void igl::opengl::ViewerCore::toggle(unsigned int &property_mask) const
{
  property_mask ^= id;
}

IGL_INLINE bool igl::opengl::ViewerCore::is_set(unsigned int property_mask) const
{
  return (property_mask & id);
}

IGL_INLINE igl::opengl::ViewerCore::ViewerCore()
{
  // Default colors
  background_color << 0.93f, 0.95f, 0.97f, 1.0f;

  // Default lights settings
  light_position_camera << 0.0f, 0.0f, 5.0f;
  //light_position << 10000.0f, 10000.0f, 10000.0f;
  light_position << 0.0f, 10000.0f, 0.0f;
  light_position_2 << -10000.0f, -10000.0f, -10000.0f;
  //light_position_3 << 0.0f, 1000.0f, 1000.0f; // unused
  //light_position_4 << 0.0f, -1000.0f, 1000.0f; // unused
  lighting_factor = 0.75f; //on

  // Default trackball
  trackball_angle = Eigen::Quaternionf::Identity();
  rotation_type = ViewerCore::ROTATION_TYPE_TRACKBALL;
  set_rotation_type(ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP);

  // Defalut model viewing parameters
  //model_zoom = 1.0f;
  //model_translation << 0, 0, 0;


  // Camera parameters
  camera_base_zoom = 1.0f;
  camera_zoom = 1.0f;
  orthographic = false;
  camera_view_angle = 45.0;
  camera_dnear = 1.0;
  camera_dfar = 100.0;
  camera_base_translation << 0, 0, 0;
  camera_translation << 0, 0, 0;
  camera_eye << 0, 0, 5;
  camera_center << 0, 0, 0;
  // camera reference frame
  camera_up << 0, 1, 0;
  camera_dir << 0, 0, -1;
  camera_side << -1, 0, 0;

  depth_test = true;

  is_animating = false;
  animation_max_fps = 60.; // was 30

  viewport.setZero();
}

IGL_INLINE void igl::opengl::ViewerCore::init()
{
}

IGL_INLINE void igl::opengl::ViewerCore::shut()
{
}
