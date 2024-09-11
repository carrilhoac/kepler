
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

void Mem::set(const char *buf, std::size_t num)
{
  dat.resize(num);
  std::memcpy(dat.data(), buf, num);
}
  
void Mem::load(const char *filename)
{
#ifdef DEBUG
  if(!std::filesystem::exists(filename))
    error("file does not exist");
#endif 

  std::size_t filesize = std::filesystem::file_size(filename);
  
#ifdef DEBUG
  if(!filesize)
    warn("empty file");
#endif 

  dat.resize(filesize);
  pos=0;
  
  std::ifstream stream(filename, std::ifstream::binary);
  stream.read(reinterpret_cast<char*>(&dat[0]), filesize);
}

void Mem::save(const char *filename) const 
{
#ifdef DEBUG
  if(!dat.size())
    warn("no data to write");
#endif 

  std::ofstream stream(filename, std::ofstream::binary);
  stream.write(reinterpret_cast<const char*>(&dat[0]), dat.size());
}

Mem& Mem::seek(std::size_t offset, eSEEKPOS from)
{
  std::size_t n=dat.size();
  
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
  dat.insert(dat.end(), b.dat.begin(), b.dat.end());
  return *this;
}

std::size_t Mem::read(char *buf, std::size_t num)
{
  std::size_t n=dat.size();
  
  if(pos+num>n)
    num=n-pos;
  
  std::memcpy(buf, &dat[0] +pos, num);
  
  pos += num;
  return num;
}

std::size_t Mem::write(const char *buf, std::size_t num)
{
  std::size_t n=dat.size();
  
  if(pos+num>n)
    num=n-pos;
  
  std::memcpy(&dat[0] + pos, buf, num);
  
  pos += num;
  return num;
}

//////////////////////////////////////////////////////////////////////
// Text utility

// length of string until control character
int strlen_ctrl(const char *s)
{
  char *p = (char*)s;
  while(!iscntrl(*p)) 
    p++;
  return p-s;
}

int count_lines(const char *s)
{
  int i,j;
  
  for(j=0,i=0;s[i];i++){
    if(s[i]=='\n'||(s[i]=='\r'&&s[i+1]!='\n')){
      j++;
    }
  }
  return j;
}

int parse_lines_n(const char *s, const char **lines, int max_lines)
{
  int i,j;
  
  lines[0]=s;  
  for(i=1,j=1;i<max_lines;j++){
    if(!s[j])
      return i;
    if(s[j]=='\n' // Win32, Unix (POSIX) and MacOS 10 or later
    ||(s[j]=='\r'&&s[j+1]!='\n')){// for MacOS 9 and older
      lines[i++]=s+(++j); // new line
    }
  }
  return i;
}

const char **parse_lines(const char *s)
{
  int num;
  const char **lines;
  
  num=count_lines(s);
  lines=(const char**)malloc(num*sizeof(char*));
  parse_lines_n(s, lines, num);
  
  return lines;
}
