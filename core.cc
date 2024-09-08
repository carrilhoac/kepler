
#include "kepler.h"
#include <filesystem>
#include <fstream>

//////////////////////////////////////////////////////////////////////
//  Debugging 

#ifdef DEBUG
void error_(const char *from, const char *msg)
{
  std::cerr<<"[ERROR] "<<from<<": "<<msg<<std::endl;
  exit(1);
}

void warn_(const char *from, const char *msg)
{
  std::cerr<<"[WARNING] "<<from<<": "<<msg<<std::endl;
}
#endif 

//////////////////////////////////////////////////////////////////////
//  Memory chunk

void Mem::load(const char *filename)
{
#ifdef DEBUG
  if(std::filesystem::exists(filename))
    error("file does not exist");
#endif 

  std::size_t filesize = std::filesystem::file_size(filename);
  
#ifdef DEBUG
  if(!filesize)
    warn("empty file");
#endif 

  data.resize(filesize);
  pos=0;
  
  std::ifstream stream(filename, std::ifstream::binary);
  stream.read(reinterpret_cast<char*>(&data[0]), filesize);
}

void Mem::save(const char *filename) const 
{
#ifdef DEBUG
  if(!data.size())
    warn("no data to write");
#endif 

  std::ofstream stream(filename, std::ofstream::binary);
  stream.write(reinterpret_cast<const char*>(&data[0]), data.size());
}

Mem& Mem::seek(std::size_t offset, eSEEKPOS from)
{
  std::size_t n=data.size();
  
  switch(from)
  {
  case BEG:
#ifdef DEBUG
    if(offset>n)
      error("position out of range");
#endif
    pos = offset;
  break;
  
  case CUR:
#ifdef DEBUG
    if(pos+offset>n)
      error("position out of range");
#endif
    pos += offset;
  break;
  
  case END:
#ifdef DEBUG
    if(offset>n)
      error("position out of range");
#endif
    pos = n-offset;
  break;
  }
  return *this;
}

Mem& Mem::append(const Mem& b)
{
  data.insert(data.end(), b.data.begin(), b.data.end());
  return *this;
}

std::size_t Mem::read(char *buf, std::size_t num)
{
  std::size_t n=data.size();
  
  if(pos+num>n)
    num=n-pos;
  
  std::memcpy(buf, &data[0] +pos, num);
  
  pos += num;
  return num;
}

std::size_t Mem::write(const char *buf, std::size_t num)
{
  std::size_t n=data.size();
  
  if(pos+num>n)
    num=n-pos;
  
  std::memcpy(&data[0] + pos, buf, num);
  
  pos += num;
  return num;
}


//////////////////////////////////////////////////////////////////////
// ASCII multi-line text 

