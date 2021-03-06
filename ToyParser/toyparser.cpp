#include <cstdlib>
#include <vector>
#include <cassert>
#include <cctype>
#include <iostream>

struct Token
{
  enum Kind {
    Eof = 0,
    LastChar = 127,
    Int,
    Name
  };
  Kind kind;
  const char* start;
  const char* stop;
  union {
    int intval;
    const char* name;
  };
};

enum struct OpCode : uint8_t {
  HALT,
  ADD,
  SUB,
  MUL,
  DIV,
  NEG,
  LIT
};

struct Intern
{
  size_t len;
  const char* str;
};


struct Context
{
  std::vector<Intern> interns;

  Token token;  // Current token.
  const char* stream = nullptr;


  std::vector<OpCode> expression;
};


void appendOpCode(Context* context, OpCode opcode)
{
  context->expression.push_back(opcode);
}

void appendLiteral(Context* context, int32_t val)
{
  context->expression.push_back(OpCode::LIT);
  auto size = context->expression.size();
  context->expression.resize(size + 4);
  *reinterpret_cast<int32_t*>(context->expression.data() + size) = val;
}

bool staticCheckExpression(Context* context, size_t maxStack)
{
  size_t stackSize = 0;

  size_t i = 0; 
  for (; i < context->expression.size(); i++) {
    switch (context->expression[i])
    {
    case OpCode::HALT:
      if (stackSize != 1) {
        std::cerr << "Unbalanced stack.\n";
        return false;
      }
      if (i + 1 != context->expression.size()) {
        std::cerr << "Premature halt-instruction.\n";
        return false;
      }
      return true;
      break;
    case OpCode::ADD:
    case OpCode::SUB:
    case OpCode::MUL:
    case OpCode::DIV:
      if (stackSize < 2) {
        std::cerr << "Stack underflow\n";
        return false;
      }
      stackSize = stackSize - 2 + 1;
      break;

    case OpCode::NEG:
      if (stackSize < 1) {
        std::cerr << "Stack underflow\n";
        return false;
      }
      stackSize = stackSize;
      break;
    case OpCode::LIT:
      if (maxStack < stackSize + 1) {
        std::cerr << "Stack overflow\n";
        return false;
      }
      i += 4;
      stackSize = stackSize + 1;
      break;
    default:
      std::cerr << "Illegal opcode " << (int)(context->expression[i]) << "\n";
      return false;
      break;
    }
  }
  std::cerr << "No halt instruction at end of expression.\n";
  return false;
}

int runExpression(Context* context, OpCode* opcodes)
{
  std::vector<int32_t> stack;
  while (*opcodes != OpCode::HALT) {
    switch (*opcodes++)
    {
    case OpCode::ADD: {
      assert(2 <= stack.size());
      auto b = stack.back(); stack.pop_back();
      auto a = stack.back(); stack.pop_back();
      stack.push_back(a + b);
      break;
    }
    case OpCode::SUB: {
      assert(2 <= stack.size());
      auto b = stack.back(); stack.pop_back();
      auto a = stack.back(); stack.pop_back();
      stack.push_back(a - b);
      break;
    }
    case OpCode::MUL: {
      assert(2 <= stack.size());
      auto b = stack.back(); stack.pop_back();
      auto a = stack.back(); stack.pop_back();
      stack.push_back(a * b);
      break;
    }
    case OpCode::DIV: {
      assert(2 <= stack.size());
      auto b = stack.back(); stack.pop_back();
      auto a = stack.back(); stack.pop_back();
      stack.push_back(a / b);
      break;
    }
    case OpCode::NEG: {
      assert(1 <= stack.size());
      stack.back() = -stack.back();
      break;
    }
    case OpCode::LIT: {
      const auto * ptr = reinterpret_cast<int32_t*>(opcodes);
      stack.push_back(*ptr);
      opcodes += sizeof(int32_t);
      break;
    }
    default:
      std::cerr << "Unrecognized opcode: " << (int)opcodes[-1] << "\n";
      return 0;
      break;
    }
  }

  if (stack.size() == 1) {
    return stack[0];
  }
  else {
    std::cerr << "Unbalanced stack, size=" << stack.size() << " at exit.\n";
    return 0;
  }
}


void testVM(Context* context)
{
  context->expression.clear();
  appendLiteral(context, 42);
  appendOpCode(context, OpCode::HALT);
  assert(staticCheckExpression(context, 1024));
  assert(runExpression(context, context->expression.data()) == 42);

  context->expression.clear();
  appendLiteral(context, 2);
  appendLiteral(context, 42);
  appendOpCode(context, OpCode::MUL);
  appendOpCode(context, OpCode::HALT);
  assert(staticCheckExpression(context, 1024));
  assert(runExpression(context, context->expression.data()) == 2 * 42);
}


const char* internString(Context* context, const char* begin, const char* end)
{
  assert(begin <= end);
  auto len = size_t(end - begin);
  for (auto & intern : context->interns) {
    if (intern.len == len && (strncmp(begin, intern.str, len) == 0)) {
      return intern.str;
    }
  }

  auto * str = new char[len + 1];
  memcpy(str, begin, len);
  str[len] = '\0';
  context->interns.emplace_back(Intern{ len, str });
  return str;
}

const char* internString(Context* context, const char* str)
{
  const auto len = strlen(str);
  return internString(context, str, str + len);
}

void stringInterningTest(Context* context)
{
  assert(std::string(internString(context, "foo")) == "foo");
  assert(internString(context, "foo") == internString(context, "foo"));
  assert(internString(context, "foo") != internString(context, "bar"));
}

void inline consumeWhitespace(Context* context)
{
  while (true)
  {
    switch (*context->stream)
    {
    case ' ':
    case '\t':
    case '\n':
    case '\v':
    case '\f':
    case '\r':
      context->stream++;
      break;
    default:
      return;
    }
  }
}

void nextToken(Context* context)
{
  auto & tok = context->token;
  auto & stream = context->stream;

  consumeWhitespace(context);
  tok.start = context->stream;
  switch (*stream)
  {
  case 'a': case 'b': case 'c': case 'd': case 'e': case 'f': case 'g': case 'h':
  case 'i': case 'j': case 'k': case 'l': case 'm': case 'n': case 'o': case 'p':
  case 'q': case 'r': case 's': case 't': case 'u': case 'v': case 'w': case 'x':
  case 'y': case 'z':
  case 'A': case 'B': case 'C': case 'D': case 'E': case 'F': case 'G': case 'H':
  case 'I': case 'J': case 'K': case 'L': case 'M': case 'N': case 'O': case 'P':
  case 'Q': case 'R': case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
  case 'Y': case 'Z': case '_':
    do {
      stream++;
    } while (std::isalnum(*stream) || *stream == '_');
    tok.kind = Token::Name;
    tok.name = internString(context, tok.start, stream);
    break;
  case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9': {
    int val = *stream++ - '0';
    while (std::isdigit(*stream)) {
      val = 10*val + *stream++ - '0';
    }
    tok.kind = Token::Int;
    tok.intval = val;
    break;
  }
  default:
    tok.kind = Token::Kind(*stream++);
    break;
  }
  tok.stop = stream;
}

void initStream(Context* context, const char* stream)
{
  context->stream = stream;
  nextToken(context);
}

bool inline isToken(Context* context, Token::Kind kind)
{
  return context->token.kind == kind;
}

bool matchToken(Context* context, Token::Kind kind)
{
  if (isToken(context, kind)) {
    nextToken(context);
    return true;
  }
  else {
    return false;
  }
}

bool expectToken(Context* context, Token::Kind kind)
{
  if (isToken(context, kind)) {
    nextToken(context);
    return true;
  }
  else {
    std::cerr << "Expected token " << (int)kind << ", got " << (int)context->token.kind << "\n";
    return false;
  }

}

#define assert_token(c,x) assert(matchToken((c), (Token::Kind)(x)))
#define assert_token_name(c,x) assert((c)->token.name == internString((c),(x)) && matchToken((c), Token::Name))
#define assert_token_int(c,x) assert((c)->token.intval == (x) && matchToken((c), Token::Int))
#define assert_token_eof(c) assert((c)->token.kind == Token::Eof);
void lexTest(Context* context)
{
  auto str = "XY + (  XY )          _HELLO1,234+994";
  initStream(context, str);
  assert_token_name(context, "XY");
  assert_token(context, '+');
  assert_token(context, '(');
  assert_token_name(context, "XY");
  assert_token(context, ')');
  assert_token_name(context, "_HELLO1");
  assert_token(context, ',');
  assert_token_int(context, 234);
  assert_token(context, '+');
  assert_token_int(context, 994);
  assert_token_eof(context);
}
#undef assert_token
#undef assert_token_name
#undef assert_token_int
#undef assert_token_eof

// expr3 = INT | '(' expr ')'
// expr2 = '-' expr2 | expr3
// expr1 = expr2([*/ ] expr2)*
// expr0 = expr1([+-] expr1)*
// expr = expr0

int parseExpr(Context* context);

int parseExpr3(Context* context)
{
  if (isToken(context, Token::Int)) {
    auto val = context->token.intval;
    appendLiteral(context, val);
    nextToken(context);
    return val;
  }
  else if (matchToken(context, Token::Kind('('))) {
    auto val = parseExpr(context);
    expectToken(context, Token::Kind(')'));
    return val;
  }
  else {
    std::cerr << "Expected integer or '(', got " << (int)context->token.kind << "\n";
    return 0;
  }

}

int parseExpr2(Context* context)
{
  if (matchToken(context, Token::Kind('-'))) {
    auto val = -parseExpr2(context);
    appendOpCode(context, OpCode::NEG);
    return val;
  }
  else if (matchToken(context, Token::Kind('+'))) {
    return parseExpr2(context);
  }
  else {
    return parseExpr3(context);
  }
}

int parseExpr1(Context* context)
{
  int val = parseExpr2(context);
  while (isToken(context, Token::Kind('*')) || isToken(context, Token::Kind('*'))) {
    auto op = context->token.kind;
    nextToken(context);
    int rval = parseExpr2(context);
    if (op == Token::Kind('*')) {
      appendOpCode(context, OpCode::MUL);
      val *= rval;
    }
    else {
      assert(op == Token::Kind('/'));
      appendOpCode(context, OpCode::DIV);
      val /= rval;
    }
  }
  return val;
}

int parseExpr0(Context* context)
{
  int val = parseExpr1(context);
  while (isToken(context, Token::Kind('+')) || isToken(context, Token::Kind('-'))) {
    auto op = context->token.kind;
    nextToken(context);
    int rval = parseExpr1(context);
    if (op == Token::Kind('+')) {
      appendOpCode(context, OpCode::ADD);
      val += rval;
    }
    else {
      assert(op == Token::Kind('-'));
      appendOpCode(context, OpCode::SUB);
      val -= rval;
    }
  }
  return val;
}

int parseExpr(Context* context)
{
  return parseExpr0(context);
}

int parseExpr(Context* context, const char* str)
{
  context->expression.clear();
  initStream(context, str);
  auto vaL = parseExpr(context);
  appendOpCode(context, OpCode::HALT);
  staticCheckExpression(context, 1024);
  return vaL;
}

#define assert_expr(c,x) assert(parseExpr((c),#x) == (x) && runExpression(context, (c)->expression.data()) == (x))
void parseTest(Context* context)
{
  assert_expr(context, 1);
  assert_expr(context, (1));
  assert_expr(context, -+1);
  assert_expr(context, 1 - 2 - 3);
  assert_expr(context, 2 * 3 + 4 * 5);
  assert_expr(context, 2 * (3 + 4) * 5);
  assert_expr(context, 2 + -3);
}
#undef assert_expr


int main(int argc, char** argv)
{
  Context context;

  stringInterningTest(&context);
  lexTest(&context);
  parseTest(&context);
  testVM(&context);

  return 0;
}

