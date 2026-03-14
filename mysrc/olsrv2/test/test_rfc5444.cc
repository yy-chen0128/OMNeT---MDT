//
// Created by GuoYudi on 2026/3/11.
//


#include <cassert>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

#include "../include/rfc5444/rfc5444.h"

using mysrc::olsrv2::rfc5444::IcvCrypt7182;
using mysrc::olsrv2::rfc5444::IcvHash7182;
using mysrc::olsrv2::rfc5444::Result;
using mysrc::olsrv2::rfc5444::rfc5444_seqno_difference;
using mysrc::olsrv2::rfc5444::rfc5497_timetlv_decode;
using mysrc::olsrv2::rfc5444::rfc5497_timetlv_encode;
using mysrc::olsrv2::rfc5444::rfc5497_timetlv_get_from_vector;
using mysrc::olsrv2::rfc5444::rfc7181_metric_decode;
using mysrc::olsrv2::rfc5444::rfc7181_metric_encode;
using mysrc::olsrv2::rfc5444::rfc7182_get_crypt_id;
using mysrc::olsrv2::rfc5444::rfc7182_get_crypt_name;
using mysrc::olsrv2::rfc5444::rfc7182_get_hash_id;
using mysrc::olsrv2::rfc5444::rfc7182_get_hash_name;
using mysrc::olsrv2::rfc5444::rfc7181_metric_field;
using mysrc::olsrv2::rfc5444::RFC7181_METRIC_MAX;
using mysrc::olsrv2::rfc5444::RFC7181_METRIC_MIN;
using mysrc::olsrv2::rfc5444::toString;

static void test_timetlv_vector() {
  // 测试函数：根据 hopcount 在 TLV 向量中查找对应的时间值。
  const std::vector<uint8_t> vec = {10, 2, 20, 5, 30};
  assert(rfc5497_timetlv_get_from_vector(vec.data(), vec.size(), 1) == 10);
  assert(rfc5497_timetlv_get_from_vector(vec.data(), vec.size(), 3) == 20);
  assert(rfc5497_timetlv_get_from_vector(vec.data(), vec.size(), 6) == 30);

  const std::vector<uint8_t> bad_vec = {1, 2};
  assert(rfc5497_timetlv_get_from_vector(bad_vec.data(), bad_vec.size(), 1) == 255);
}

static void test_timetlv_encode_decode() {
  // 测试函数：RFC5497 时间值的编码与解码。
  assert(rfc5497_timetlv_encode(0) == 0);
  assert(rfc5497_timetlv_decode(0) == 0);

  const uint64_t input_ms = 2000;
  const uint8_t encoded = rfc5497_timetlv_encode(input_ms);
  const uint64_t decoded = rfc5497_timetlv_decode(encoded);
  assert(decoded >= 1000);
}

static void test_metric_encode_decode() {
  // 测试函数：RFC7181 度量值的编码与解码。
  rfc7181_metric_field field{};
  const uint32_t metric = RFC7181_METRIC_MIN + 42;
  assert(rfc7181_metric_encode(&field, metric) == 0);
  const uint32_t decoded = rfc7181_metric_decode(&field);
  assert(decoded == metric);

  rfc7181_metric_field high_field{};
  assert(rfc7181_metric_encode(&high_field, RFC7181_METRIC_MAX + 1) == 0);
  assert(high_field.b[0] == 0x0f);
  assert(high_field.b[1] == 0xff);
}

static void test_seqno_difference() {
  // 测试函数：计算 RFC5444 序列号差值（含回绕）。
  assert(rfc5444_seqno_difference(10, 9) == 1);
  assert(rfc5444_seqno_difference(1, 65535) == 2);
}

static void test_iana_helpers() {
  // 测试函数：RFC7182 哈希/加密名称与 ID 的互相转换。
  assert(std::string_view(rfc7182_get_hash_name(IcvHash7182::Sha256)) == "sha256");
  assert(rfc7182_get_hash_id("SHA256") == IcvHash7182::Sha256);
  assert(std::string_view(rfc7182_get_crypt_name(IcvCrypt7182::Aes)) == "aes");
  assert(rfc7182_get_crypt_id("AES") == IcvCrypt7182::Aes);
}

static void test_result_to_string() {
  // 测试函数：结果码到字符串的映射。
  assert(std::string_view(toString(Result::Okay)) == "Okay");
  assert(std::string_view(toString(Result::BadTlvIndex)) == "BadTlvIndex");
}

int main() {
  test_timetlv_vector();
  test_timetlv_encode_decode();
  test_metric_encode_decode();
  test_seqno_difference();
  test_iana_helpers();
  test_result_to_string();

  std::cout << "rfc5444 tests passed\n";
  return 0;
}
