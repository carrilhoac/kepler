
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
// ASCII multi-line text 


Text::Text(const char *s)
{
  dat.set(s, strlen(s));
  parse();
}

Text& Text::operator=(const char *s)
{
  dat.set(s, strlen(s));
  parse();
  return *this;
}
  
const char* Text::line(std::size_t index) const
{
#ifdef DEBUG
  if(index>=lines.size())
    error("position out of range");
#endif
  return lines[index];
}

char* Text::line(std::size_t index)
{
#ifdef DEBUG
  if(index>=lines.size())
    error("position out of range");
#endif
  return lines[index];
}

void Text::load(const char *filepath)
{
  dat.load(filepath);
  parse();
}

void Text::save(const char *filepath, const char *eol) const
{
  std::size_t n=lines.size();
  std::ofstream stream(filepath, std::ofstream::binary);
  
  for(std::size_t i=0; i<n; i++)
    stream<< lines[i] << eol;
}
  
int Text::find_line(const char *key, int start_line) const
{
  std::size_t n=lines.size();
  std::size_t m=strlen(key);
  
  for(std::size_t i=start_line; i<n; i++)
    if(strlen(lines[i])>=m)
      if(strstr(lines[i],key))
        return i;
  
  return -1;
}

std::ostream& operator<<(std::ostream& os, const Text& t)
{
  std::size_t n=t.lines.size();
  
  for(std::size_t i=0; i<n; i++)
    os<<t.lines[i]<<std::endl;
  
  return os;
}

void Text::parse()
{
  std::size_t n=dat.size();
  char *buf=dat.data();
    
  lines.clear();
  
  if(!n)
    return;
  
  lines.push_back(&buf[0]);
  for(std::size_t i=1; i<n-1; i++){
    if(buf[i]=='\n' || (buf[i]=='\r'&&buf[i+1]!='\n')){
      buf[i]=0;
      lines.push_back(&buf[i+1]);
    }
  }
}
  