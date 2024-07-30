<template>
  <view class="container">
    <button @click="startRecording">开始录音</button>
    <button @click="stopRecording">停止录音</button>
	{{pitchData}}
    <view v-if="pitchData.length > 0">
      <text v-for="(pitch, index) in pitchData" :key="index">音高: {{ pitch }} Hz</text>
    </view>
  </view>
</template>

<script>
import PitchFinder from 'pitchfinder';
const detectPitch = PitchFinder.AMDF();

export default {
  data() {
    return {
      recorderManager: null,
      pitchData: [],
    };
  },
  methods: {
    startRecording() {
      this.recorderManager = uni.getRecorderManager();
	  
	console.log('oooo', this.recorderManager)
	  
      this.recorderManager.onFrameRecorded((res) => {
        const pcmData = new Float32Array(res.frameBuffer);
		console.error('处理音频：', pcmData)
        const pitch = detectPitch(pcmData);
        if (pitch) {
          this.pitchData.push(pitch);
        }
      });
	console.error('12222', this.recorderManager.start)
      this.recorderManager.start();
	  console.error(('2222'))
    },
    stopRecording() {
      if (this.recorderManager) {
        this.recorderManager.stop();
      }
    },
  },
};
</script>

<style>
.container {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100vh;
}
button {
  margin: 10px;
}
</style>
