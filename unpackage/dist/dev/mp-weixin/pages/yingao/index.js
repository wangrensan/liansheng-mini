"use strict";
const common_vendor = require("../../common/vendor.js");
const detectPitch = common_vendor._default.AMDF();
const _sfc_main = {
  data() {
    return {
      recorderManager: null,
      pitchData: []
    };
  },
  methods: {
    startRecording() {
      this.recorderManager = common_vendor.index.getRecorderManager();
      console.log("oooo", this.recorderManager);
      this.recorderManager.onFrameRecorded((res) => {
        const pcmData = new Float32Array(res.frameBuffer);
        console.error("处理音频：", pcmData);
        const pitch = detectPitch(pcmData);
        if (pitch) {
          this.pitchData.push(pitch);
        }
      });
      console.error("12222", this.recorderManager.start);
      this.recorderManager.start();
      console.error("2222");
    },
    stopRecording() {
      if (this.recorderManager) {
        this.recorderManager.stop();
      }
    }
  }
};
function _sfc_render(_ctx, _cache, $props, $setup, $data, $options) {
  return common_vendor.e({
    a: common_vendor.o((...args) => $options.startRecording && $options.startRecording(...args)),
    b: common_vendor.o((...args) => $options.stopRecording && $options.stopRecording(...args)),
    c: common_vendor.t($data.pitchData),
    d: $data.pitchData.length > 0
  }, $data.pitchData.length > 0 ? {
    e: common_vendor.f($data.pitchData, (pitch, index, i0) => {
      return {
        a: common_vendor.t(pitch),
        b: index
      };
    })
  } : {});
}
const MiniProgramPage = /* @__PURE__ */ common_vendor._export_sfc(_sfc_main, [["render", _sfc_render]]);
wx.createPage(MiniProgramPage);
