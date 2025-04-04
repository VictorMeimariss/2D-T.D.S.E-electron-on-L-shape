import { AutoModel, AutoProcessor, RawImage } from '@xenova/transformers';

// Get URL from command-line argument
const url = process.argv[2];

if (!url)
{
    console.error("Error, URL was not provided.");
    process.exit(1);
}

console.log("Processing URL:", url);

// Load model
const model = await AutoModel.from_pretrained('Xenova/gelan-c', {
    // quantized: false,    // (Optional) Use unquantized version.
})

// Load processor
const processor = await AutoProcessor.from_pretrained('Xenova/gelan-c');
// processor.feature_extractor.do_resize = false;                   // (Optional) Disable resizing
// processor.feature_extractor.size = { width: 128, height: 128 }   // (Optional) Update resize value

// Read image and run processor
const image = await RawImage.read(url);
const { pixel_values } = await processor(image);

// Run object detection
const { outputs } = await model({ images: pixel_values })
const predictions = outputs.tolist();

for (const [xmin, ymin, xmax, ymax, score, id] of predictions) {
    const bbox = [xmin, ymin, xmax, ymax].map(x => x.toFixed(2)).join(', ')
    console.log(`Found "${model.config.id2label[id]}" at [${bbox}] with score ${score.toFixed(2)}.`)
}
