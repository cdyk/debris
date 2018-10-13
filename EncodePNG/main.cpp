#if 0
// FIXME: add zlib dependency.
#include <vector>
#include <zlib.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <time.h>
#include <xmmintrin.h>
#include <smmintrin.h>

#define WIDTH 1024 //6
#define HEIGHT 768

class TimeStamp
{
public:
  TimeStamp()
  {
    clock_gettime(CLOCK_MONOTONIC, &m_time);
  }

  static
    double
    delta(const TimeStamp& start, const TimeStamp& stop)
  {
    timespec delta;
    if ((stop.m_time.tv_nsec - start.m_time.tv_nsec) < 0) {
      delta.tv_sec = stop.m_time.tv_sec - start.m_time.tv_sec - 1;
      delta.tv_nsec = 1000000000 + stop.m_time.tv_nsec - start.m_time.tv_nsec;
    }
    else {
      delta.tv_sec = stop.m_time.tv_sec - start.m_time.tv_sec;
      delta.tv_nsec = stop.m_time.tv_nsec - start.m_time.tv_nsec;
    }
    return delta.tv_sec + 1e-9*delta.tv_nsec;
  }


protected:
  timespec m_time;
};




class BitPusher
{
public:
  BitPusher(std::vector<unsigned char>& data)
    : m_data(data),
    m_pending_data(0u),
    m_pending_count(0u)
  {}

  ~BitPusher()
  {
    while (m_pending_count >= 8u) {
      m_data.push_back(m_pending_data);
      m_pending_data = m_pending_data >> 8u;
      m_pending_count -= 8u;
    }
    if (m_pending_count) {
      m_data.push_back(m_pending_data);
    }
  }

  void
    pushBits(unsigned long long int bits, unsigned int count)
  {
    m_pending_data = m_pending_data | (bits << m_pending_count);
    m_pending_count += count;
    while (m_pending_count >= 8u) {
      m_data.push_back(m_pending_data);
      m_pending_data = m_pending_data >> 8u;
      m_pending_count -= 8u;
    }
  }

  void
    pushBitsReverse(unsigned long long int bits, unsigned int count)
  {
    //std::cerr << "pushing " << count << " bits:\t";
    //for(int i=0; i<count; i++) {
    //     std::cerr << ((bits>>(count-i-1))&1);
    //}
    //std::cerr << "\n";

    bits = bits << (64 - count);
    bits = ((bits & 0x5555555555555555ull) << 1) | ((bits >> 1) & 0x5555555555555555ull);
    bits = ((bits & 0x3333333333333333ull) << 2) | ((bits >> 2) & 0x3333333333333333ull);
    bits = ((bits & 0x0F0F0F0F0F0F0F0Full) << 4) | ((bits >> 4) & 0x0F0F0F0F0F0F0F0Full);
    bits = ((bits & 0x00FF00FF00FF00FFull) << 8) | ((bits >> 8) & 0x00FF00FF00FF00FFull);
    bits = ((bits & 0x0000FFFF0000FFFFull) << 16) | ((bits >> 16) & 0x0000FFFF0000FFFFull);
    bits = ((bits & 0x00000000FFFFFFFFull) << 32) | ((bits >> 32) & 0x00000000FFFFFFFFull);
    pushBits(bits, count);
  }


protected:
  std::vector<unsigned char>& m_data;
  unsigned long long int      m_pending_data;
  unsigned int                m_pending_count;


};


void
encodeCount(unsigned int& bits, unsigned int& bits_n, unsigned int count)
{
  //std::cerr << "count=" << count << "  ";


  uint q = 0;

  if (count < 3) {
    abort();
    // no copies, 1 or 2 never occurs.
  }
  else if (count < 11) {
    unsigned int t = count - 3;
    t = (t & 0x07u) + ((257 - 256) << 0);
    bits = (bits << 7) | t;           // 7-bit count code
    bits_n += 7;
  }
  else if (count < 19) {     // count=11..18, code=265..268, 1 extra bit
    int t = count - 11;
    q = (t & 0x01u) << (8 - 1);
    t = (t & 0x06u) + ((265 - 256) << 1);
    bits = (bits << 8) | t;
    bits_n += 8;
  }
  else if (count < 35) {     // count=19..34, code=269..272, 2 extra bits
    int t = count - 19;
    q = (t & 0x03u) << (8 - 2);
    t = (t & 0x0Cu) + ((269 - 256) << 2);
    bits = (bits << 9) | t;
    bits_n += 9;
  }
  else if (count < 67) {     // c=35..66, code=273..276, 3 extra bits
    int t = count - 35;
    q = (t & 0x07u) << (8 - 3);
    t = (t & 0x18u) + ((273 - 256) << 3);
    bits = (bits << 10) | t;
    bits_n += 10;
  }
  else if (count < 115) {    // c=67..114, 7-bit code=277..279, 4 extra bits
    int t = count - 67;
    q = (t & 0x0Fu) << (8 - 4);
    t = (t & 0x30u) + ((277 - 256) << 4);
    bits = (bits << 11) | t;
    bits_n += 11;
  }
  else if (count < 131) {    // c=115..130, 8-bit code=280, 4 extra bits
    int t = count - 115;
    q = (t & 0x0Fu) << (8 - 4);
    t = (t & 0x30u) + ((280 - 280 + 0xc0) << 4);
    bits = (bits << 12) | t;
    bits_n += 12;
  }
  else if (count < 258) {    // c=131..257, code=281..284, 5 extra bits
    int t = count - 131;
    bits = (bits << 12) | ((t&(0xff << 5)) + ((281 - 280 + 0xc0) << 5));
    q = (t & 0x1Fu) << (8 - 5);
    bits_n += 13;
  }
  else if (count < 259) {
    bits = (bits << 8) | (285 - 280 + 0xc0);
    bits_n += 8;
  }
  else {
    std::cerr << "unsupported count " << count << "\n";
    abort();
  }

  q = ((q & 0x55u) << 1) | ((q >> 1) & 0x55u);
  q = ((q & 0x33u) << 2) | ((q >> 2) & 0x33u);
  q = ((q & 0x0Fu) << 4) | ((q >> 4) & 0x0Fu);
  bits = bits | q;
}


void
encodeDistance(unsigned int& bits, unsigned int& bits_n, unsigned int distance)
{
  //std::cerr << "distance=" << distance << "  ";

  uint q = 0;

  if (distance < 1) {
    std::cerr << "unsupported distance " << distance << "\n";
    abort();
  }
  else if (distance < 5) {   // d=1..4, code=0..3, 0 extra bits
    bits = (bits << 5) | (distance - 1);
    bits_n += 5;
  }
  else if (distance < 9) {   // d=5..8, code=4..5, 1 extra bit
    unsigned int t = distance - 5;
    bits = (bits << 6) | ((t&(0x1f << 1)) + (4 << 1));
    q = (t & 1u) << (32 - 1);
    bits_n += 6;
  }
  else if (distance < 17) {      // d=9..16, code=6..7, 2 extra bits
    unsigned int t = distance - 9;
    bits = (bits << 7) | ((t&(0x1 << 2)) + (6 << 2));
    q = (t & 3u) << (32 - 2);
    bits_n += 7;
  }
  else if (distance < 33) {      // d=17..32, code=8..9, 3 extra bits
    unsigned int t = distance - 17;
    bits = (bits << 8) | ((t&(0x1 << 3)) + (8 << 3));
    q = (t & 7u) << (32 - 3);
    bits_n += 8;
  }
  else if (distance < 65) {      // d=33..64, code=10..11, 4 extra bits
    unsigned int t = distance - 33;
    bits = (bits << 9) | ((t&(0x1 << 4)) + (10 << 4));
    q = (t & 0xFu) << (32 - 4);
    bits_n += 9;
  }
  else if (distance < 129) {     // d=65..128, code=12..13, 5 extra bits
    unsigned int t = distance - 65;
    bits = (bits << 10) | ((t&(0x1 << 5)) + (12 << 5));
    q = (t & 0x1Fu) << (32 - 5);
    bits_n += 10;
  }
  else if (distance < 257) {     // d=129..256, code=14,15, 6 extra bits
    unsigned int t = distance - 129;
    q = (t & 0x3Fu) << (32 - 6);
    t = (t & 0x40u) + (14 << 6);
    bits = (bits << 11) | t;
    bits_n += 11;
  }
  else if (distance < 513) {     // d=257..512, code 16..17, 7 extra bits
    unsigned int t = distance - 257;
    q = (t & 0x7Fu) << (32 - 7);
    t = (t & 0x80u) + (16 << 7);
    bits = (bits << 12) | t;
    bits_n += 12;
  }
  else if (distance < 1025) {     // d=257..512, code 18..19, 8 extra bits
    unsigned int t = distance - 513;
    q = (t & 0x0FFu) << (32 - 8);
    t = (t & 0x100u) + (18 << 8);
    bits = (bits << 13) | t;
    bits_n += 13;
  }
  else if (distance < 2049) {    // d=1025..2048, code 20..21, 9 extra bits
    unsigned int t = distance - 1025;
    q = (t & 0x1FFu) << (32 - 9);
    t = (t & 0x200u) + (20 << 9);
    bits = (bits << 14) | t;
    bits_n += 14;
  }
  else if (distance < 4097) {    // d=2049..4096, code 22..23, 10 extra bits
    unsigned int t = distance - 2049;
    q = (t & 0x3FFu) << (32 - 10);
    t = (t & 0x400u) + (22 << 10);
    bits = (bits << 15) | t;
    bits_n += 15;
  }
  else if (distance < 8193) {    // d=4097..8192, code 24..25, 11 extra bits
    unsigned int t = distance - 4097;
    q = (t & 0x7FFu) << (32 - 11);
    t = (t & 0x800u) + (24 << 11);
    bits = (bits << 16) | t;
    bits_n += 16;
  }
  else if (distance < 16385) {   // d=8193..16384, code 26..27, 12 extra bits
    unsigned int t = distance - 8193;
    q = (t & 0x0FFFu) << (32 - 12);
    t = (t & 0x1000u) + (26 << 12);
    bits = (bits << 17) | t;
    bits_n += 17;
  }
  else if (distance < 32769) {   // d=16385..32768, code 28..29, 13 extra bits
    unsigned int t = distance - 16385;
    q = (t & 0x1FFFu) << (32 - 13);
    t = (t & 0x2000u) + (28 << 13);
    bits = (bits << 18) | t;
    bits_n += 18;
  }
  else {
    std::cerr << "Illegal " << distance << "\n";
    abort();
  }

  q = ((q & 0x55555555u) << 1) | ((q >> 1) & 0x55555555u);
  q = ((q & 0x33333333u) << 2) | ((q >> 2) & 0x33333333u);
  q = ((q & 0x0F0F0F0Fu) << 4) | ((q >> 4) & 0x0F0F0F0Fu);
  q = ((q & 0x00FF00FFu) << 8) | ((q >> 8) & 0x00FF00FFu);
  q = ((q & 0x0000FFFFu) << 16) | ((q >> 16) & 0x0000FFFFu);
  bits = bits | q;
}




void
encodeLiteralTriplet(unsigned int& bits, unsigned int& bits_n, unsigned int rgb)
{
  //std::cerr << "literal=" << rgb << "  ";

  unsigned int t = (rgb >> 16) & 0xffu;
  if (t < 144) {
    bits = (t + 48);
    bits_n = 8;
  }
  else {
    bits = (t + 256);
    bits_n = 9;
  }

  t = (rgb >> 8) & 0xffu;
  if (t < 144) {
    bits = (bits << 8) | (t + 48);
    bits_n += 8;
  }
  else {
    bits = (bits << 9) | (t + 256);
    bits_n += 9;
  }

  t = (rgb) & 0xffu;
  if (t < 144) {
    bits = (bits << 8) | (t + 48);
    bits_n += 8;
  }
  else {
    bits = (bits << 9) | (t + 256);
    bits_n += 9;
  }
}


unsigned long int
CRC(const std::vector<unsigned long>& crc_table, const unsigned char* p, size_t length)
{
  size_t i;
  unsigned long int crc = 0xffffffffl;
  for (i = 0; i < length; i++) {
    unsigned int ix = (p[i] ^ crc) & 0xff;
    crc = crc_table[ix] ^ ((crc >> 8));
  }
  return ~crc;
}




void
writeIDAT3(std::ofstream& file, const std::vector<unsigned char>& img, const std::vector<unsigned long>& crc_table)
{
  TimeStamp start;


  std::vector<unsigned char> IDAT(8);
  // IDAT chunk header
  IDAT[4] = 'I';
  IDAT[5] = 'D';
  IDAT[6] = 'A';
  IDAT[7] = 'T';

  // --- create deflate chunk ------------------------------------------------
  IDAT.push_back(8 + (7 << 4));  // CM=8=deflate, CINFO=7=32K window size = 112
  IDAT.push_back(94 /* 28*/);           // FLG

  unsigned int dat_size;


  unsigned int s1 = 1;
  unsigned int s2 = 0;
  {
    BitPusher pusher(IDAT);
    pusher.pushBitsReverse(6, 3);    // 5 = 101

    std::vector<unsigned int> buffer(WIDTH * 5, ~0u);

    unsigned int* zrows[4] = {
        buffer.data(),
        buffer.data() + WIDTH,
        buffer.data() + 2 * WIDTH,
        buffer.data() + 3 * WIDTH,
    };

    unsigned int o = 0;


    for (int j = 0; j < HEIGHT; j++) {
      unsigned int* rows[4] = {
          zrows[(j + 3) & 1],
          zrows[(j + 2) & 1],
          zrows[(j + 1) & 1],
          zrows[j & 1]
      };


      //            unsigned int* p_row = rows[ j&1 ];
      //unsigned int* c_row = rows[ (j+1)&1 ];

      pusher.pushBitsReverse(0 + 48, 8);    // Push png scanline filter
      s1 += 0;                                // update adler 1 & 2
      s2 += s1;

      int match_src_i = 0;
      int match_src_j = 0;
      int match_dst_i = 0;
      int match_dst_o = 0; // used for debugging
      int match_length = 0;
      //std::cerr << o << "---\n";
      o++;

      for (int i = 0; i < WIDTH; i++) {
        unsigned int R = img[3 * (WIDTH*j + i) + 0];
        unsigned int G = img[3 * (WIDTH*j + i) + 1];
        unsigned int B = img[3 * (WIDTH*j + i) + 2];
        unsigned int RGB = (R << 16) | (G << 8) | B;

        s1 = (s1 + R);
        s2 = (s2 + s1);
        s1 = (s1 + G);
        s2 = (s2 + s1);
        s1 = (s1 + B);
        s2 = (s2 + s1);

        //std::cerr << o << ":\t" << RGB << "\t";
        rows[0][i] = RGB;

        bool redo;

        do {
          redo = false;
          bool emit_match = false;
          bool emit_verbatim = false;

          if (match_length == 0) {
            // We have no running matching, try to find candidate
            int k;


            match_src_j = 0;
            for (k = i - 1; (k >= 0) && (rows[0][k] != RGB); k--) {}

            if (k < 0) {
              match_src_j = 1;
              for (k = WIDTH - 1; (k >= 0) && (rows[1][k] != RGB); k--) {}
            }

            if (k >= 0) {
              // Found match,
              match_src_i = k;
              match_dst_i = i;
              match_dst_o = o;
              match_length = 1;
              if (i == WIDTH - 1) {
                emit_match = true;
              }
            }
            else {
              emit_verbatim = true;
            }
          }
          else {
            // We are matching

            if (match_length >= 86) {
              // Max matchlength is 86*3=258, flush and continue
              emit_match = true;
              redo = true;
            }

            else if ((match_src_i + match_length < WIDTH) &&    // don't match outside scanline
              (rows[match_src_j][match_src_i + match_length] == RGB))  // check if matching
            {
              match_length = match_length + 1;
              if (i == WIDTH - 1) {
                emit_match = true;
              }
            }
            else {
              emit_match = true;
              redo = true;

              // try to find new match source
              int k = match_src_i - 1;
              for (int m = match_src_j; (emit_match) && (m < 2); m++) {
                for (; (emit_match) && (k >= 0); k--) {
                  bool fail = false;
                  for (int l = 0; l <= match_length; l++) {
                    if (rows[0][match_dst_i + l] != rows[m][k + l]) {
                      fail = true;
                      break;
                    }
                  }
                  if (!fail) {
                    match_src_j = m;
                    match_src_i = k;
                    match_length = match_length + 1;
                    emit_match = false;
                    redo = false;
                    break;
                  }
                }
                k = WIDTH - 1;
              }

            }
          }

          if (((i == WIDTH - 1) && (match_length)) || emit_match) {
            unsigned int bits = 0;
            unsigned int bits_n = 0;

            unsigned int count = 3 * match_length;
            encodeCount(bits, bits_n, count);
            pusher.pushBitsReverse(bits, bits_n);

            bits = 0;
            bits_n = 0;
            unsigned int distance = 3 * (match_dst_i - match_src_i);

            if (match_src_j > 0) {
              distance += 3 * WIDTH + 1;
            }

            encodeDistance(bits, bits_n, distance);
            pusher.pushBitsReverse(bits, bits_n);
            match_length = 0;
          }

          if (emit_verbatim) {
            unsigned int bits = 0;
            unsigned int bits_n = 0;
            encodeLiteralTriplet(bits, bits_n, RGB);
            pusher.pushBitsReverse(bits, bits_n);
          }
        } while (redo);

        o += 3;
      }
      s1 = s1 % 65521;
      s2 = s2 % 65521;

    }
    pusher.pushBits(0, 7);    // EOB 
  }
  unsigned int adler = (s2 << 16) + s1;

  IDAT.push_back(((adler) >> 24) & 0xffu); // Adler32
  IDAT.push_back(((adler) >> 16) & 0xffu);
  IDAT.push_back(((adler) >> 8) & 0xffu);
  IDAT.push_back(((adler) >> 0) & 0xffu);

  // --- end deflate chunk --------------------------------------------------


  // Update PNG chunk content size for IDAT
  dat_size = IDAT.size() - 8u;
  IDAT[0] = ((dat_size) >> 24) & 0xffu;
  IDAT[1] = ((dat_size) >> 16) & 0xffu;
  IDAT[2] = ((dat_size) >> 8) & 0xffu;
  IDAT[3] = ((dat_size) >> 0) & 0xffu;

  unsigned long crc = CRC(crc_table, IDAT.data() + 4, dat_size + 4);
  IDAT.resize(IDAT.size() + 4u);  // make room for CRC
  IDAT[dat_size + 8] = ((crc) >> 24) & 0xffu;
  IDAT[dat_size + 9] = ((crc) >> 16) & 0xffu;
  IDAT[dat_size + 10] = ((crc) >> 8) & 0xffu;
  IDAT[dat_size + 11] = ((crc) >> 0) & 0xffu;


  TimeStamp stop;
  std::cerr << "IDAT3 used " << TimeStamp::delta(start, stop) << "\n";


  if (1) {
    std::vector<unsigned char> quux(10 * 1024 * 1024);

    z_stream stream;
    int err;

    stream.next_in = (z_const Bytef *)IDAT.data() + 8;
    stream.avail_in = (uInt)IDAT.size() - 8;
    stream.next_out = quux.data();
    stream.avail_out = quux.size();
    stream.zalloc = (alloc_func)0;
    stream.zfree = (free_func)0;

    err = inflateInit(&stream);
    if (err != Z_OK) {
      std::cerr << "inflateInit failed: " << err << "\n";
      abort();
    }

    err = inflate(&stream, Z_FINISH);

    if (stream.msg != NULL) {
      std::cerr << stream.msg << "\n";
    }


    uLongf quux_size = quux.size();
    err = uncompress(quux.data(), &quux_size, IDAT.data() + 8, IDAT.size() - 8);
    if (err != Z_OK) {
      std::cerr << "uncompress="
        << err
        << "\n";
    }
    if (quux_size != ((3 * WIDTH + 1)*HEIGHT)) {
      std::cerr << "uncompress_size=" << quux_size
        << ", should be=" << ((3 * WIDTH + 1)*HEIGHT) << "\n";
    }
  }


  file.write(reinterpret_cast<char*>(IDAT.data()), dat_size + 12);
}





void
createCRCTable(std::vector<unsigned long>& crc_table)
{
  crc_table.resize(256);
  for (int j = 0; j < 256; j++) {
    unsigned long int c = j;
    for (int i = 0; i < 8; i++) {
      if (c & 0x1) {
        c = 0xedb88320ul ^ (c >> 1);
      }
      else {
        c = c >> 1;
      }
    }
    crc_table[j] = c;
  }
}

void
createDummyImage(std::vector<unsigned char>& img)
{
  img.resize(3 * WIDTH*HEIGHT);
  for (int j = 0; j < HEIGHT; j++) {
    for (int i = 0; i < WIDTH; i++) {
      float x = i / (WIDTH / 2.f) - 1.f;
      float y = j / (HEIGHT / 2.f) - 1.f;
      float r = sqrt(x*x + y * y);
      if (r < 0.9f) {
        img[3 * (WIDTH*j + i) + 0] = 255;
        img[3 * (WIDTH*j + i) + 1] = 255 * 0.5f*(sinf(30 * r) + 1.f);;
        img[3 * (WIDTH*j + i) + 2] = 255 * 0.5f*(sinf(40 * r) + 1.f);
      }
      else if (r < 1.f) {
        int q = (int)(200 * r);
        if (q & 1) {
          img[3 * (WIDTH*j + i) + 0] = q;
          img[3 * (WIDTH*j + i) + 1] = 255;
          img[3 * (WIDTH*j + i) + 2] = 255;
        }
        else {
          img[3 * (WIDTH*j + i) + 0] = 0;
          img[3 * (WIDTH*j + i) + 1] = q;
          img[3 * (WIDTH*j + i) + 2] = 0;
        }
      }
      else {
        img[3 * (WIDTH*j + i) + 0] = 255;
        img[3 * (WIDTH*j + i) + 1] = 255;
        img[3 * (WIDTH*j + i) + 2] = 0;
      }
    }
  }
}

void
writeSignature(std::ofstream& file)
{
  unsigned char signature[8] =
  {
      137, 80, 78, 71, 13, 10, 26, 10
  };
  file.write(reinterpret_cast<char*>(signature), sizeof(signature));
}


void
writeIHDR(std::ofstream& file, const std::vector<unsigned long>& crc_table)
{
  // IHDR chunk, 13 + 12 (length, type, crc) = 25 bytes
  unsigned char IHDR[25] =
  {
    // Chunk length (4 bytes)
    ((13) >> 24) & 0xffu,((13) >> 16) & 0xffu, ((13) >> 8) & 0xffu, ((13) >> 0) & 0xffu,
    // Chunk type (4 bytes)
    'I', 'H', 'D', 'R',
    // Image width (4 bytes)
    ((WIDTH) >> 24) & 0xffu,((WIDTH) >> 16) & 0xffu, ((WIDTH) >> 8) & 0xffu, ((WIDTH) >> 0) & 0xffu,
    // Image height
    ((HEIGHT) >> 24) & 0xffu,((HEIGHT) >> 16) & 0xffu, ((HEIGHT) >> 8) & 0xffu, ((HEIGHT) >> 0) & 0xffu,
    // bits per channel, RGB triple, ..., .., image not interlaced (5 bytes)
    8, 2, 0, 0, 0,
    // CRC of 13+4 bytes
    0, 0, 0, 0
  };
  unsigned long crc = CRC(crc_table, IHDR + 4, 13 + 4);
  IHDR[21] = ((crc) >> 24) & 0xffu;    // image width
  IHDR[22] = ((crc) >> 16) & 0xffu;
  IHDR[23] = ((crc) >> 8) & 0xffu;
  IHDR[24] = ((crc) >> 0) & 0xffu;

  file.write(reinterpret_cast<char*>(IHDR), sizeof(IHDR));
}


void
writeIDAT(std::ofstream& file, const std::vector<unsigned char>& img, const std::vector<unsigned long>& crc_table)
{
  TimeStamp start;

  std::vector<unsigned char> filtered((3 * WIDTH + 1)*HEIGHT);
  for (int j = 0; j < HEIGHT; j++) {
    filtered[(3 * WIDTH + 1)*j + 0] = 0;
    for (int i = 0; i < WIDTH; i++) {
      filtered[(3 * WIDTH + 1)*j + 1 + 3 * i + 0] = img[3 * WIDTH*j + 3 * i + 0];
      filtered[(3 * WIDTH + 1)*j + 1 + 3 * i + 1] = img[3 * WIDTH*j + 3 * i + 1];
      filtered[(3 * WIDTH + 1)*j + 1 + 3 * i + 2] = img[3 * WIDTH*j + 3 * i + 2];
    }
  }
  std::vector<unsigned char> IDAT(16 * 1024 * 1024);

  uLongf dat_size = IDAT.size() - 12;

  TimeStamp c_start;
  int c = compress((Bytef*)(IDAT.data() + 8), &dat_size, (Bytef*)filtered.data(), filtered.size());
  TimeStamp c_stop;


  if (c == Z_MEM_ERROR) {
    std::cerr << "Z_MEM_ERROR\n";
    exit(EXIT_FAILURE);
  }
  else if (c == Z_BUF_ERROR) {
    std::cerr << "Z_BUF_ERROR\n";
    exit(EXIT_FAILURE);
  }

  IDAT[0] = ((dat_size) >> 24) & 0xffu;
  IDAT[1] = ((dat_size) >> 16) & 0xffu;
  IDAT[2] = ((dat_size) >> 8) & 0xffu;
  IDAT[3] = ((dat_size) >> 0) & 0xffu;
  IDAT[4] = 'I';
  IDAT[5] = 'D';
  IDAT[6] = 'A';
  IDAT[7] = 'T';

  unsigned long crc = CRC(crc_table, IDAT.data() + 4, dat_size + 4);
  IDAT[dat_size + 8] = ((crc) >> 24) & 0xffu;    // image width
  IDAT[dat_size + 9] = ((crc) >> 16) & 0xffu;
  IDAT[dat_size + 10] = ((crc) >> 8) & 0xffu;
  IDAT[dat_size + 11] = ((crc) >> 0) & 0xffu;

  TimeStamp stop;
  std::cerr << "IDAT used " << TimeStamp::delta(start, stop)
    << ", compress used " << TimeStamp::delta(c_start, c_stop) << "\n";

  file.write(reinterpret_cast<char*>(IDAT.data()), dat_size + 12);
}


void
writeIDAT2(std::ofstream& file, const std::vector<unsigned char>& img, const std::vector<unsigned long>& crc_table)
{
  TimeStamp start;


  std::vector<unsigned char> IDAT(8);
  // IDAT chunk header
  IDAT[4] = 'I';
  IDAT[5] = 'D';
  IDAT[6] = 'A';
  IDAT[7] = 'T';

  // --- create deflate chunk ------------------------------------------------
  IDAT.push_back(8 + (7 << 4));  // CM=8=deflate, CINFO=7=32K window size = 112
  IDAT.push_back(94 /* 28*/);           // FLG

  unsigned int dat_size;
  unsigned int s1 = 1;
  unsigned int s2 = 0;
  {
    BitPusher pusher(IDAT);
    pusher.pushBitsReverse(6, 3);    // 5 = 101

    for (int j = 0; j < HEIGHT; j++) {

      // push scan-line filter type
      pusher.pushBitsReverse(1 + 48, 8);    // Push a 1 (diff with left)
      s1 += 1;                                // update adler 1 & 2
      s2 += s1;

#if 1
      unsigned int trgb_p = 0xffffffff;
      unsigned int c = 0;
      for (int i = 0; i < WIDTH; i++) {

        unsigned int trgb_l = 0;
        for (int k = 0; k < 3; k++) {
          unsigned int t = (img[3 * WIDTH*j + 3 * i + k]
                            - (i == 0 ? 0 : img[3 * WIDTH*j + 3 * (i - 1) + k])) & 0xffu;

          s1 = (s1 + t);                  // update adler 1 & 2
          s2 = (s2 + s1);
          trgb_l = (trgb_l << 8) | t;
        }
        if ((i == 0) || (i == WIDTH - 1) || (trgb_l != trgb_p) || (c >= 66)) {
          // flush copies
          if (c == 0) {
            // no copies, 1 or 2 never occurs.
          }
          else if (c < 11) {
            pusher.pushBitsReverse(c - 2, 7);
            pusher.pushBitsReverse(2, 5);
          }
          else if (c < 19) {
            int t = c - 11;
            pusher.pushBitsReverse((t >> 1) + 9, 7);
            pusher.pushBitsReverse((t & 1), 1);
            pusher.pushBitsReverse(2, 5);
          }
          else if (c < 35) {
            int t = c - 19;
            pusher.pushBitsReverse((t >> 2) + 13, 7);
            pusher.pushBits((t & 3), 2);
            pusher.pushBitsReverse(2, 5);
          }
          else if (c < 67) {
            int t = c - 35;
            pusher.pushBitsReverse((t >> 3) + 17, 7);
            pusher.pushBits((t & 7), 3);
            pusher.pushBitsReverse(2, 5);
          }
          c = 0;


          // need to write literal
          int r = (trgb_l >> 16) & 0xffu;
          if (r < 144) {
            pusher.pushBitsReverse(r + 48, 8);
          }
          else {
            pusher.pushBitsReverse(r + (400 - 144), 9);
          }
          int g = (trgb_l >> 8) & 0xffu;
          if (g < 144) {
            pusher.pushBitsReverse(g + 48, 8);
          }
          else {
            pusher.pushBitsReverse(g + (400 - 144), 9);
          }
          int b = (trgb_l >> 0) & 0xffu;
          if (b < 144) {
            pusher.pushBitsReverse(b + 48, 8);
          }
          else {
            pusher.pushBitsReverse(b + (400 - 144), 9);
          }
          trgb_p = trgb_l;
        }
        else {
          c += 3;
        }
      }
#else            
      for (int i = 0; i < WIDTH; i++) {





        for (int k = 0; k < 3; k++) {
          int q = img[3 * WIDTH*j + 3 * i + k];
          s1 = (s1 + q);                  // update adler 1 & 2
          s2 = (s2 + s1);

          if (q < 144) {
            pusher.pushBitsReverse(q + 48, 8);            // 48 = 00110000
          }
          else {
            pusher.pushBitsReverse(q + (400 - 144), 9);     // 400 = 110010000
          }
        }
      }
#endif

      // We can do up to 5552 iterations before we need to run the modulo,
      // see comment on NMAX adler32.c in zlib.
      s1 = s1 % 65521;
      s2 = s2 % 65521;
    }
    pusher.pushBits(0, 7);    // EOB 
  }
  unsigned int adler = (s2 << 16) + s1;

  IDAT.push_back(((adler) >> 24) & 0xffu); // Adler32
  IDAT.push_back(((adler) >> 16) & 0xffu);
  IDAT.push_back(((adler) >> 8) & 0xffu);
  IDAT.push_back(((adler) >> 0) & 0xffu);

  // --- end deflate chunk --------------------------------------------------




  // Update PNG chunk content size for IDAT
  dat_size = IDAT.size() - 8u;
  IDAT[0] = ((dat_size) >> 24) & 0xffu;
  IDAT[1] = ((dat_size) >> 16) & 0xffu;
  IDAT[2] = ((dat_size) >> 8) & 0xffu;
  IDAT[3] = ((dat_size) >> 0) & 0xffu;

  unsigned long crc = CRC(crc_table, IDAT.data() + 4, dat_size + 4);
  IDAT.resize(IDAT.size() + 4u);  // make room for CRC
  IDAT[dat_size + 8] = ((crc) >> 24) & 0xffu;
  IDAT[dat_size + 9] = ((crc) >> 16) & 0xffu;
  IDAT[dat_size + 10] = ((crc) >> 8) & 0xffu;
  IDAT[dat_size + 11] = ((crc) >> 0) & 0xffu;


  TimeStamp stop;
  std::cerr << "IDAT2 used " << TimeStamp::delta(start, stop) << "\n";


  if (1) {
    std::vector<unsigned char> quux(10 * 1024 * 1024);

    z_stream stream;
    int err;

    stream.next_in = (z_const Bytef *)IDAT.data() + 8;
    stream.avail_in = (uInt)IDAT.size() - 8;
    stream.next_out = quux.data();
    stream.avail_out = quux.size();
    stream.zalloc = (alloc_func)0;
    stream.zfree = (free_func)0;

    err = inflateInit(&stream);
    if (err != Z_OK) {
      std::cerr << "inflateInit failed: " << err << "\n";
      abort();
    }

    err = inflate(&stream, Z_FINISH);

    if (stream.msg != NULL) {
      std::cerr << stream.msg << "\n";
    }


    uLongf quux_size = quux.size();
    err = uncompress(quux.data(), &quux_size, IDAT.data() + 8, IDAT.size() - 8);
    if (err != Z_OK) {
      std::cerr << "uncompress="
        << err
        << "\n";
    }
    if (quux_size != ((3 * WIDTH + 1)*HEIGHT)) {
      std::cerr << "uncompress_size=" << quux_size
        << ", should be=" << ((3 * WIDTH + 1)*HEIGHT) << "\n";
    }

  }


  file.write(reinterpret_cast<char*>(IDAT.data()), dat_size + 12);
}



void
writeIEND(std::ofstream& file, const std::vector<unsigned long>& crc_table)
{
  unsigned char IEND[12] = {
      0, 0, 0, 0,         // payload size
      'I', 'E', 'N', 'D', // chunk id
      174, 66, 96, 130    // chunk crc
  };
  file.write(reinterpret_cast<char*>(IEND), 12);
}

#endif
int
main(int argc, char **argv)
{
#if 0
  std::vector<unsigned long> crc_table;
  createCRCTable(crc_table);

  std::vector<unsigned char> img;
  createDummyImage(img);

  std::ofstream png("img.png");
  writeSignature(png);
  writeIHDR(png, crc_table);
  writeIDAT(png, img, crc_table);
  writeIEND(png, crc_table);
  png.close();

  std::ofstream png2("img2.png");
  writeSignature(png2);
  writeIHDR(png2, crc_table);
  writeIDAT2(png2, img, crc_table);
  writeIEND(png2, crc_table);
  png2.close();

  std::ofstream png3("img3.png");
  writeSignature(png3);
  writeIHDR(png3, crc_table);
  writeIDAT3(png3, img, crc_table);
  writeIEND(png3, crc_table);
  png3.close();

  {
    std::ofstream png2("img2.png");
    writeSignature(png2);
    writeIHDR(png2, crc_table);
    writeIDAT2(png2, img, crc_table);
    writeIEND(png2, crc_table);
    png2.close();
  }
#endif

  return 0;
}
