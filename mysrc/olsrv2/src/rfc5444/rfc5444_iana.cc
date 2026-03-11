#include <array>
#include <cstddef>
#include <string_view>
#include <strings.h>

#include "../../include/rfc5444/rfc5444.h"

namespace mysrc::olsrv2::rfc5444 {
namespace {
using namespace std::literals;

constexpr std::string_view kUnknown = "unknown"sv;

constexpr std::size_t kHashCount = static_cast<std::size_t>(IcvHash7182::Count);
constexpr std::size_t kCryptCount = static_cast<std::size_t>(IcvCrypt7182::Count);

constexpr std::array<std::string_view, kHashCount> kHashes = [] {
  std::array<std::string_view, kHashCount> a{};
  a[static_cast<std::size_t>(IcvHash7182::Identity)] = "identity"sv;
  a[static_cast<std::size_t>(IcvHash7182::Sha1)] = "sha1"sv;
  a[static_cast<std::size_t>(IcvHash7182::Sha224)] = "sha224"sv;
  a[static_cast<std::size_t>(IcvHash7182::Sha256)] = "sha256"sv;
  a[static_cast<std::size_t>(IcvHash7182::Sha384)] = "sha384"sv;
  a[static_cast<std::size_t>(IcvHash7182::Sha512)] = "sha512"sv;
  return a;
}();

constexpr std::array<std::string_view, kCryptCount> kCrypt = [] {
  std::array<std::string_view, kCryptCount> a{};
  a[static_cast<std::size_t>(IcvCrypt7182::Identity)] = "identity"sv;
  a[static_cast<std::size_t>(IcvCrypt7182::Rsa)] = "rsa"sv;
  a[static_cast<std::size_t>(IcvCrypt7182::Dsa)] = "dsa"sv;
  a[static_cast<std::size_t>(IcvCrypt7182::Hmac)] = "hmac"sv;
  a[static_cast<std::size_t>(IcvCrypt7182::Des3)] = "3des"sv;
  a[static_cast<std::size_t>(IcvCrypt7182::Aes)] = "aes"sv;
  a[static_cast<std::size_t>(IcvCrypt7182::Ecdsa)] = "ecdsa"sv;
  return a;
}();

// Bridge arrays for legacy API: `const char **`
std::array<const char *, kHashCount> gHashCstr{};
std::array<const char *, kCryptCount> gCryptCstr{};

inline const char *sv_to_cstr(std::string_view sv) {
  return sv.empty() ? kUnknown.data() : sv.data();
}

inline void ensure_bridges_initialized() {
  static bool initialized = false;
  if (initialized) return;

  for (std::size_t i = 0; i < gHashCstr.size(); ++i) {
    gHashCstr[i] = sv_to_cstr(kHashes[i]);
  }
  for (std::size_t i = 0; i < gCryptCstr.size(); ++i) {
    gCryptCstr[i] = sv_to_cstr(kCrypt[i]);
  }
  initialized = true;
}
} // namespace

const char *rfc7182_get_hash_name(IcvHash7182 hash) {
  const int idx = static_cast<int>(hash);
  if (idx < 0 || idx >= static_cast<int>(kHashes.size())) {
    return kUnknown.data();
  }
  return sv_to_cstr(kHashes[static_cast<std::size_t>(idx)]);
}

const char *rfc7182_get_crypt_name(IcvCrypt7182 crypt) {
  const int idx = static_cast<int>(crypt);
  if (idx < 0 || idx >= static_cast<int>(kCrypt.size())) {
    return kUnknown.data();
  }
  return sv_to_cstr(kCrypt[static_cast<std::size_t>(idx)]);
}

const char **rfc7182_get_hashes(void) {
  ensure_bridges_initialized();
  return gHashCstr.data();
}

const char **rfc7182_get_crypto(void) {
  ensure_bridges_initialized();
  return gCryptCstr.data();
}

IcvHash7182 rfc7182_get_hash_id(const char *name) {
  if (name == nullptr) {
    return IcvHash7182::Unknown;
  }
  for (std::size_t i = 0; i < kHashes.size(); ++i) {
    if (kHashes[i].empty()) continue;
    if (strcasecmp(kHashes[i].data(), name) == 0) {
      return static_cast<IcvHash7182>(i);
    }
  }
  return IcvHash7182::Unknown;
}

IcvCrypt7182 rfc7182_get_crypt_id(const char *name) {
  if (name == nullptr) {
    return IcvCrypt7182::Unknown;
  }
  for (std::size_t i = 0; i < kCrypt.size(); ++i) {
    if (kCrypt[i].empty()) continue;
    if (strcasecmp(kCrypt[i].data(), name) == 0) {
      return static_cast<IcvCrypt7182>(i);
    }
  }
  return IcvCrypt7182::Unknown;
}

} // namespace mysrc::olsrv2::rfc5444

