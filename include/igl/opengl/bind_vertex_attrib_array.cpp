#include "bind_vertex_attrib_array.h"
#include <iostream>

IGL_INLINE GLint igl::opengl::bind_vertex_attrib_array(
  const GLuint program_shader,
  const std::string &name, 
  GLuint bufferID, 
  const Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &M, 
  bool refresh)
{
  GLint id = glGetAttribLocation(program_shader, name.c_str());
  //std::cout << "Line 12: id = " << id << std::endl;
  if (id < 0)
    return id;
  if (M.size() == 0)
  {
    glDisableVertexAttribArray(id);
    return id;
  }
  glBindBuffer(GL_ARRAY_BUFFER, bufferID);
  
  if (refresh) {
      //std::cout << "de data staat nu in the videokaart" << std::endl;
      glBufferData(GL_ARRAY_BUFFER, sizeof(float) * M.size(), M.data(), GL_DYNAMIC_DRAW);
  }
  glVertexAttribPointer(id, M.cols(), GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(id);
  return id;
}

IGL_INLINE GLint igl::opengl::bind_vertex_attrib_array_print_init(
    const GLuint program_shader,
    const std::string& name,
    GLuint bufferID,
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M,
    bool refresh,
    const GLuint n_rows)
{    
    GLint id = glGetAttribLocation(program_shader, name.c_str());
    if (id < 0)
        return id;
    if (M.size() == 0)
    {
        //std::cout << "Er is geen data voor " << name << std::endl;
        glDisableVertexAttribArray(id);
        return id;
    }
    glBindBuffer(GL_ARRAY_BUFFER, bufferID);
    if (refresh) {
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * M.size() * n_rows, nullptr, GL_DYNAMIC_DRAW); // initialize buffer
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * M.size(), M.data()); // set first data block in buffer 
        //std::cout << "de data van " << name << " staat nu in the videokaart" << std::endl;
    }
    glVertexAttribPointer(id, M.cols(), GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(id);
    return id;
}

IGL_INLINE GLint igl::opengl::bind_vertex_attrib_array_print(
    const GLuint program_shader,
    const std::string& name,
    GLuint bufferID,
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M,
    bool refresh,
    GLuint step)
{
    GLint id = glGetAttribLocation(program_shader, name.c_str());
    if (id < 0)
        return id;
    if (M.size() == 0)
    {
        //std::cout << "Er is geen data voor " << name << std::endl;
        glDisableVertexAttribArray(id);
        return id;
    }
    glBindBuffer(GL_ARRAY_BUFFER, bufferID);
    //if (refresh)
        //glBufferData(GL_ARRAY_BUFFER, sizeof(float) * M.size(), M.data(), GL_DYNAMIC_DRAW);
    if (step >= 5)
        if (refresh) {
            glBufferSubData(GL_ARRAY_BUFFER, (sizeof(float) * M.size()*(step-4)), sizeof(float) * M.size(), M.data()); // set data block in buffer 
            //std::cout << "de data van " << name << " staat nu in the videokaart" << std::endl;
        }
            
    glVertexAttribPointer(id, M.cols(), GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(id);
    return id;
}
